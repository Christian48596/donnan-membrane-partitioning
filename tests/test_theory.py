"""Independent checks of conservation, equilibrium, activity scales and rank."""
import sys
from pathlib import Path
import unittest
import numpy as np
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from hpm.manning import STOICH, gamma_counterion, gamma_coion, solve_mobile_salt
from hpm.pitzer import (SALTS, gamma_pm, gamma_pm_molar, molality_from_molarity,
                        molar_activity_product, mean_activity_product, WATER_DENSITY_KG_L)
from hpm.properties import solution_gamma, solution_property
from hpm.thermodynamics import Species, solve_donnan, mobile_ionic_strength, association_constant
from hpm.association import (external_species, external_association, paired_partition,
                             manning_species, CHARGES)
from hpm.exchange import (population_weight, separation_factor, ceiling, required_K,
                          required_energy, solve_finite_bath, condensed_fraction)
from hpm.baseline_error import Measurement, check_identifiability, simulate_power, design_matrix, _rank
from hpm.moments import moment1, moment2, exact_ln_K
from hpm.screening import debye_length, crossover_concentration
import analysis
import dle


class ManningEquilibrium(unittest.TestCase):
    def test_primary_denominator_and_coefficients(self):
        # Independent numerical values from the equations, not copied closure code.
        additions={"NaCl":2.,"MgCl2":6.,"Na2SO4":6.,"MgSO4":8.}
        for salt,n in additions.items():
            s=STOICH[salt]; X,xi=4.,3.1; A=X/(2*(X+n*xi))
            self.assertAlmostEqual(gamma_coion(X,s,xi),np.exp(-A*(s.z_c/s.z_g)**2),14)
            pre=(X/(s.z_g*xi)+s.z_g*s.nu_g)/(X+s.z_g*s.nu_g)
            self.assertAlmostEqual(gamma_counterion(X,s,xi),pre*np.exp(-A),14)

    def test_general_neutral_salt_product_and_charge(self):
        for salt in ("NaCl","MgCl2","Na2SO4","MgSO4"):
            s=STOICH[salt]
            for CA in (0.,3.7):
                for C in (.01,.3,1.):
                    xi=3.1; g=solution_gamma(salt,C)
                    Cm=solve_mobile_salt(salt,CA,C,g,xi,h_g=1.3,h_c=.7)
                    cg,cc=CA/s.z_g+s.nu_g*Cm,s.nu_c*Cm
                    gg,gc=gamma_counterion(CA/Cm,s,xi),gamma_coion(CA/Cm,s,xi)
                    self.assertAlmostEqual(s.z_g*cg-s.z_c*cc,CA,11)
                    lhs=(gg*cg)**s.nu_g*(gc*cc)**s.nu_c
                    rhs=1.3**s.nu_g*.7**s.nu_c*g**(s.nu_g+s.nu_c)*(s.nu_g*C)**s.nu_g*(s.nu_c*C)**s.nu_c
                    self.assertLess(abs(lhs/rhs-1),2e-10)
                    # Both single-ion relations must imply the SAME Donnan potential.
                    psi_g=-np.log(gg*cg/(1.3*g*s.nu_g*C))/s.z_g
                    psi_c=np.log(gc*cc/(.7*g*s.nu_c*C))/s.z_c
                    self.assertAlmostEqual(psi_g,psi_c,10)

    def test_activity_limits_and_condensation_domain(self):
        for s in STOICH.values():
            self.assertEqual(gamma_coion(0,s,3.1),1.)
            self.assertEqual(gamma_counterion(0,s,3.1),1.)
            self.assertAlmostEqual(gamma_counterion(1e12,s,3.1),np.exp(-.5)/(s.z_g*3.1),10)
            self.assertAlmostEqual(gamma_coion(1e12,s,3.1),np.exp(-.5*(s.z_c/s.z_g)**2),9)
        with self.assertRaises(ValueError): gamma_coion(1,STOICH['NaCl'],.9)
        self.assertEqual(solve_mobile_salt("NaCl",4.,0.,1.,3.1),0.)


class ActivityAndCharge(unittest.TestCase):
    def test_molal_to_molar_activity_identity(self):
        for salt in ("NaCl","MgCl2","Na2SO4","MgSO4"):
            for c in (.01,.1,.3,1.):
                rho=solution_property(salt,c); m=molality_from_molarity(salt,c,rho)
                gc=gamma_pm_molar(salt,c,rho)
                self.assertAlmostEqual(gc*c,WATER_DENSITY_KG_L*gamma_pm(salt,m)*m,13)
                nu=SALTS[salt].nu
                self.assertAlmostEqual(molar_activity_product(salt,c,rho)/
                                       (mean_activity_product(salt,m)*WATER_DENSITY_KG_L**nu),1.,12)

    def test_negative_concentrations_rejected(self):
        with self.assertRaises(ValueError): gamma_pm("NaCl",-.1)
        with self.assertRaises(ValueError): molality_from_molarity("NaCl",1.,.01)
        with self.assertRaises(ValueError): solution_property("MgCl2",1.,strict=True)

    def test_species_donnan_charge_and_neutral_transfer(self):
        cs={"Na+":.2,"SO4-2":.08,"NaSO4-":.04,"neutral":.1}
        sp={n:Species(n,z,h=1.4 if n=='neutral' else 1.) for n,z in
            (("Na+",1),("SO4-2",-2),("NaSO4-",-1),("neutral",0))}
        psi,cm=solve_donnan(cs,sp,3.5)
        self.assertLess(abs(sum(sp[n].z*c for n,c in cm.items())-3.5),1e-10)
        self.assertAlmostEqual(cm['neutral']/cs['neutral'],1.4)
        for n in cs: self.assertAlmostEqual(np.log(cm[n]/cs[n]),sp[n].log_prefactor-sp[n].z*psi,11)

    def test_material_ionic_strength(self):
        ii=analysis.interior_ionic_strength()
        for r in ii.itertuples():
            s=STOICH[r.salt]
            co=analysis.sorption_data().query('membrane==@r.membrane and salt_system==@r.salt and c_external_M==@r.c_ext').iloc[0].c_coion_M
            self.assertAlmostEqual(r.I_mobile_dissociated_M,.5*(s.z_g**2*r.counterion_total_M+s.z_c**2*co),12)
            self.assertLessEqual(r.I_mobile_min_association_M,r.I_mobile_dissociated_M)
        self.assertAlmostEqual(mobile_ionic_strength({'Na+':4.,'Cl-':.1,'MgSO4':9.},
                                                    {'Na+':1,'Cl-':-1,'MgSO4':0}),2.05)

    def test_apparent_charge_recovers_known_ideal_and_flags_opposite_sign(self):
        for salt in ("NaCl","MgCl2","Na2SO4","MgSO4"):
            co=analysis.ideal_donnan_coion(salt,.1,4.)
            self.assertAlmostEqual(dle.fit_fixed_charge(salt,.1,co),4.,9)
        ca=dle.apparent_charge_table()
        self.assertTrue((ca.C_apparent<0).any())
        self.assertTrue((ca.admissible_negative_membrane==(ca.C_apparent>=0)).all())

    def test_crossover_is_debye_inverse(self):
        I=crossover_concentration(78.3)
        self.assertAlmostEqual(I,1.025615,5)
        self.assertAlmostEqual(float(debye_length(I,78.3))/0.30e-9,1.,12)


class SulfateSpecies(unittest.TestCase):
    def test_external_element_and_charge_balances(self):
        for salt in ("Na2SO4","MgSO4"):
            for a in (0.,.3,.8):
                c=external_species(salt,.4,a)
                self.assertAlmostEqual(sum(CHARGES[n]*v for n,v in c.items()),0.,14)
                if salt=='Na2SO4': self.assertAlmostEqual(c['Na+']+c['NaSO4-'],.8)
                else: self.assertAlmostEqual(c['Mg2+']+c['MgSO4'],.4)
                P='NaSO4-' if salt=='Na2SO4' else 'MgSO4'
                self.assertAlmostEqual(c['SO4-2']+c[P],.4)

    def test_distinct_mass_action_stoichiometry(self):
        ideal=lambda c:{n:1. for n in c}
        for K in (.1,5.,170.):
            c=.12; t=K*c
            aMg=1-2/(1+np.sqrt(1+4*t))
            # Stable small root of t a^2 -(3t+1)a+2t=0.
            aNa=4*t/(3*t+1+np.sqrt((3*t+1)**2-8*t*t))
            self.assertAlmostEqual(external_association('MgSO4',c,K,ideal)[0],aMg,11)
            self.assertAlmostEqual(external_association('Na2SO4',c,K,ideal)[0],aNa,11)

    def test_binary_limit_of_general_manning(self):
        for salt,M,P in (("Na2SO4","Na+","NaSO4-"),("MgSO4","Mg2+","MgSO4")):
            s=STOICH[salt]; CA,C,xi=3.7,.09,3.1
            cm={M:CA/s.z_g+s.nu_g*C,'SO4-2':C,P:0.}
            gm=manning_species(cm,CA,xi,M)
            self.assertAlmostEqual(gm[M],gamma_counterion(CA/C,s,xi),12)
            self.assertAlmostEqual(gm['SO4-2'],gamma_coion(CA/C,s,xi),12)

    def test_species_equilibrium_charge_and_association_cycles(self):
        for salt in ('Na2SO4','MgSO4'):
            for c in (.01,.1,1.):
                for xi in (1.83,3.1):
                    for hp in (.7,1.,1.4):
                        p=paired_partition(salt,c,3.8,xi,h_pair=hp)
                        self.assertLess(abs(p['charge_error_M']),1e-11)
                        self.assertLess(abs(p['cycle_error']),1e-9)
                        for n,cm in p['membrane'].items():
                            cs=p['solution'][n]; gm=p['gamma_membrane'][n]; gs=p['gamma_solution'][n]
                            h=hp if n in ('MgSO4','NaSO4-') else 1.
                            self.assertAlmostEqual(np.log(gm*cm/(gs*cs)),np.log(h)-CHARGES[n]*p['psi_D'],8)

    def test_neutrality_does_not_impose_unit_partition(self):
        p=paired_partition('MgSO4',.1,3.7,3.1,h_pair=1.4,gamma_pair_membrane=2.)
        self.assertAlmostEqual(p['membrane']['MgSO4']/p['solution']['MgSO4'],.7,12)
        self.assertLess(abs(p['cycle_error']),1e-9)


class LithiumExchange(unittest.TestCase):
    def test_unselective_and_free_preference_limits(self):
        for xi in (.8,1.83,3.1):
            for y in (.01,.1,.5,.9):
                self.assertAlmostEqual(separation_factor(1.,xi,y=y),1.,12)
                self.assertAlmostEqual(separation_factor(1.,xi,y=y,q=2.),2.,11)

    def test_feed_and_salt_dependent_asymptote(self):
        for y in (.01,.1,.5):
            for co in (0.,.5,2.):
                w=population_weight(3.1,1.,co)
                expected=1+w/(y*(1-w))
                self.assertAlmostEqual(ceiling(3.1,y=y,coion_M=co),expected,10)
                self.assertLess(abs(separation_factor(1e14,3.1,y=y,coion_M=co)/expected-1),1e-10)
        self.assertAlmostEqual(ceiling(3.1),5.2)
        self.assertAlmostEqual(ceiling(1.83),2.66)
        self.assertLess(ceiling(3.1,coion_M=.5),ceiling(3.1))

    def test_exchange_inversion_and_energy(self):
        for y in (.01,.1,.5):
            for target in (1.5,2.,3.,5.):
                K=required_K(target,3.1,y=y)
                self.assertAlmostEqual(separation_factor(K,3.1,y=y),target,10)
        for a,dg in ((2.,-2.671),(3.,-4.692),(5.,-11.969)):
            self.assertAlmostEqual(required_energy(a,3.1),dg,3)
        self.assertIsNone(required_K(10.,3.1))

    def test_finite_bath_element_conservation_and_lithium_inventory(self):
        r=solve_finite_bath(.02,1.38,1.,.1,4.,0.,3.1,1e8)
        self.assertAlmostEqual(r['Li_bath_mol']+r['Li_membrane_mol'],.02,11)
        self.assertAlmostEqual(r['Na_bath_mol']+r['Na_membrane_mol'],1.38,11)
        self.assertLess(r['x'],.02/(condensed_fraction(3.1)*4.*.1)+1e-10)
        self.assertLess(r['y_final'],.02)


class IdentifiabilityAndMoments(unittest.TestCase):
    def test_fixed_counterion_valence_variance_is_confoundable(self):
        p=(Measurement('NaCl','Cl-',-1,.181,1),Measurement('Na2SO4','SO4-2',-2,.230,1),
           Measurement('LiCl','anionB',-1,.25,1))
        self.assertEqual(check_identifiability(p)['variance_added_rank'],0)
        self.assertTrue(np.isnan(simulate_power(p)['power']))

    def test_available_coions_do_not_gain_an_independent_variance_direction(self):
        for salts in (('NaCl','MgCl2','Na2SO4'),('NaCl','MgCl2','Na2SO4','MgSO4')):
            p=analysis._panel(salts)
            self.assertEqual(check_identifiability(p)['variance_added_rank'],0)
            self.assertEqual(_rank(design_matrix(p,5,False,False)),_rank(design_matrix(p,5,True,False)))

    def test_prospective_full_panel_and_false_positive_rate(self):
        p=analysis._panel(('NaCl','MgCl2','Na2SO4','MgSO4'),True)
        self.assertEqual(check_identifiability(p)['variance_added_rank'],4)
        null=simulate_power(p,amplitude=0.,mu_amplitude=.5,n_trials=3000,seed=24)
        self.assertGreater(null['power'],.03);self.assertLess(null['power'],.07)

    def test_closed_moments_and_exact_partition_match_independent_quadrature(self):
        for z in (-2,1,2):
            for k0 in (0.,2*np.pi/.65e-9):
                rp,lam,psi,theta=1.e-9,.5e-9,2.,.3
                profile=lambda t:psi*np.exp(-rp*t/lam)*np.cos(k0*rp*t+theta)
                M1=quad(lambda t:2*(1-t)*profile(t),0,1,epsabs=1e-12)[0]
                M2=quad(lambda t:2*(1-t)*profile(t)**2,0,1,epsabs=1e-12)[0]
                logK=np.log(quad(lambda t:2*(1-t)*np.exp(-z*profile(t)),0,1,epsabs=1e-12)[0])
                self.assertAlmostEqual(moment1(psi,lam,k0,rp,theta),M1,11)
                self.assertAlmostEqual(moment2(psi,lam,k0,rp,theta),M2,11)
                self.assertAlmostEqual(exact_ln_K(z,psi,lam,k0,rp,theta),logK,11)
                self.assertGreaterEqual(logK+z*M1,-1e-11) # Jensen's inequality

if __name__=='__main__': unittest.main()
