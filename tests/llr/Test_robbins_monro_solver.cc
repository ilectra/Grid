/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_llr_robbinsmonro.cc

Copyright (C) 2015, 2026

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: Guido Cossu <guido.cossu@ed.ac.uk>
Author: Ryan Hill <guido.cossu@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>
#include <Grid/qcd/llr/RobbinsMonroSolverModule.h>

int main(int argc, char **argv) 
{
  using namespace Grid;

  // Ignore command line arguments and enforce --grid 4.4.4.4
  int updatedArgc = 3;
  char** updatedArgv = (char**)malloc(sizeof(char*) * 4);
  updatedArgv[0] = argv[0];
  updatedArgv[1] = (char*)(new std::string("--grid"))->c_str();
  updatedArgv[2] = (char*)(new std::string("4.4.4.4"))->c_str();
  updatedArgv[3] = nullptr;
  Grid_init(&updatedArgc, &updatedArgv);
  GridLogLayout();

  // Starting and expected final values for solver parameters
  RealD alpha0 = 1, S0  = 1000, sigma  = 3.0;
  int n_thermalisation = 100, n_rethermalisation = 10, n_accumulation = 5, n_iter = 14;
  int n_traj = n_thermalisation + n_accumulation + (n_iter -1) * (n_accumulation + n_rethermalisation); // 300 for default values
  RealD a_expected = alpha0;
  RealD epsilon_a = 0.3;
  int n_steps = 8;
  std::vector<RealD> a_init_scan(n_steps);
  std::string gaugeGroup = "SU(3)";
  if (Sp2n_config) {
    a_expected = 9;
    a_init_scan = {-10, -1, 0, 1, 7, 9, 11, 20};
    gaugeGroup = "Sp(4)";
  }
  else {
    a_expected = 4.55;
    a_init_scan = {-10, -1, 0, 1, 4, 4.5, 6, 20};
  }
  std::cout << GridLogMessage << "Test_robbins_monro_solver: Comparing with expected values for " << gaugeGroup << ": S0 = " << S0 << ", a = " << a_expected << std::endl;

  // Scan starting values for a (including negative and extremes) and run test
  std::vector<bool> converged(n_steps);
  for (int i=0; i<n_steps; i++) {
    // The HMC runner
    typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
    HMCWrapper TheHMC;
    
    // The grid
    TheHMC.Resources.AddFourDimGrid("gauge");
    
    // Checkpointer definition
    CheckpointerParameters CPparams;  
    CPparams.config_prefix = "ckpoint_lat";
    CPparams.rng_prefix = "ckpoint_rng";
    CPparams.saveInterval = 65535;
    CPparams.format = "IEEE64BIG";
    
    TheHMC.Resources.LoadNerscCheckpointer(CPparams);

    // The RNG
    RNGModuleParameters RNGpar;
    RNGpar.serial_seeds = "1 2 3 4 5";
    RNGpar.parallel_seeds = "6 7 8 9 10";
    TheHMC.Resources.SetRNGSeeds(RNGpar);
    
    // Construct observables
    typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
    TheHMC.Resources.AddObservable<PlaqObs>();
    
    typedef TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;
    TopologyObsParameters TopParams;
    TopParams.interval = 1;
    TopParams.do_smearing = false;
    TopParams.Smearing.init_step_size = 0.01;
    TopParams.Smearing.tolerance = 1e-5;
    //  TopParams.Smearing.steps = 200;
    //  TopParams.Smearing.step_size = 0.01;
    TopParams.Smearing.meas_interval = 50;
    TopParams.Smearing.maxTau = 2.0; 
    TheHMC.Resources.AddObservable<QObs>(TopParams);
    
    // Robbins-Monro updater
    WilsonGaugeActionR bare_action(1.0);
    RealD alpha = a_init_scan[i];
    typedef ConstrainedAction<WilsonGaugeActionR> ConstrainedWilsonGaugeAction;
    ConstrainedActionParameters action_parameters{.a=alpha, .S0=S0, .sigma=sigma};
    ConstrainedWilsonGaugeAction constrained_action(bare_action, action_parameters);
    typedef RobbinsMonroSolver<ConstrainedWilsonGaugeAction, &RM> Solver;
    Solver solver(constrained_action, RobbinsMonroParameters(n_thermalisation, n_rethermalisation, n_accumulation, 1.0));
    typedef RobbinsMonroSolverModule<Solver> RmMod;
    TheHMC.Resources.AddObservable<RmMod>(solver);
    
    // Collect actions
    ActionLevel<HMCWrapper::Field> Level1(1);
    Level1.push_back(&constrained_action);
    TheHMC.TheAction.push_back(Level1);
    
    // HMC parameters are serialisable
    TheHMC.Parameters.Trajectories = n_traj;
    TheHMC.Parameters.NoMetropolisUntil = 0;
    TheHMC.Parameters.MD.MDsteps = 100;
    TheHMC.Parameters.MD.trajL   = 1.0;
    
    //TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
    
    TheHMC.Run();  // no smearing
    
    // Report intermediate action values
    auto theObservables = TheHMC.Resources.GetObservables();
    int nObs = theObservables.size();
    typedef LLRActionLogger<Solver> ActionLogger;
    // The LLRActionLogger is the second-to-last observable and not the last one,
    // because the checkpointer is always added to the end of the observables list
    // by GetObservables().
    // This is fragile, it would be nice to be able to fish out the observable name instead
    ActionLogger* theRMmodule = reinterpret_cast<ActionLogger*>(theObservables[nObs-2]);
    std::vector<std::pair<int, RealD> > theActions = theRMmodule->getMeanActionVector();
    std::cout << GridLogMessage << "Test_robbins_monro_solver: Pairs of trajectory-mean_action for every iteration: " << theActions << std::endl;
    
    // test final values
    RealD mean_action = solver.status().last_update.mean_action;
    bool S0_converged = std::abs(mean_action - S0) < 2 * sigma;
    bool a_converged  = std::abs(constrained_action.parameters().a - a_expected) < epsilon_a;
    converged[i] = S0_converged && a_converged;
    std::cout << GridLogMessage << "Test_robbins_monro_solver: a_initial = " << alpha << (S0_converged ? ", S0 converged" : ", S0 did not converge") << (a_converged ? ", a converged" : ", a did not converge") << std::endl;
    std::cout << GridLogMessage << "Test_robbins_monro_solver: a_final = " << constrained_action.parameters().a << ", Sunconstrained = " << mean_action << std::endl;
  }
  
  assert(std::all_of(converged.begin(), converged.end(), [](bool v) { return v; }));

  Grid_finalize();

} // main
