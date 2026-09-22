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

/*! @brief Example HMC run of the LLR workflow.
 *
 * Uses the constrained action
 *   \f$\beta S[U] \rightarrow aS[U] + (S[U]-S_0)^2/(2\sigma^2)\f$
 * The initial \f$a\f$ and \f$S_0\f$ can be set on the command line with
 * `--a_init <value> --S0_init <value>`.
 */

int main(int argc, char **argv) 
{
  using namespace Grid;

  Grid_init(&argc, &argv);
  GridLogLayout();

  // Starting values for solver parameters
  RealD alpha0 = 1;
  RealD S0  = 1000;
  RealD sigma  = 3.0;
  std::string arg;
  if( GridCmdOptionExists(argv, argv+argc, "--a_init") ) {
    arg = GridCmdOptionPayload(argv, argv+argc, "--a_init");
    GridCmdOptionFloat(arg, alpha0);
  }
  if( GridCmdOptionExists(argv, argv+argc, "--S0_init") ) {
    arg = GridCmdOptionPayload(argv, argv+argc, "--S0_init");
    GridCmdOptionFloat(arg, S0);
  }
  std::string gaugeGroup = "SU(3)";
  if (Sp2n_config) {
    gaugeGroup = "Sp(4)";
  }
  std::cout << GridLogMessage << "Example_LLR: Robbins-Monro solver parameters for " << gaugeGroup << ": S0 = " << S0 << ", a = " << alpha0 << ", sigma = " << sigma << std::endl;

  // The HMC runner
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  HMCWrapper TheHMC;
  
  // Grid from the command line
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
  typedef ConstrainedAction<WilsonGaugeActionR> ConstrainedWilsonGaugeAction;
  ConstrainedActionParameters action_parameters{.a=alpha0, .S0=S0, .sigma=sigma};
  ConstrainedWilsonGaugeAction constrained_action(bare_action, action_parameters);
  typedef RobbinsMonroSolver<ConstrainedWilsonGaugeAction> Solver;
  Solver solver(constrained_action, RobbinsMonroParameters(100, 10, 5, 1.0));
  typedef RobbinsMonroSolverModule<Solver> RmMod;
  TheHMC.Resources.AddObservable<RmMod>(solver);
    
  // Collect actions
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&constrained_action);
  TheHMC.TheAction.push_back(Level1);
    
  // HMC parameters are serialisable
  TheHMC.Parameters.Trajectories = 305;
  TheHMC.Parameters.NoMetropolisUntil = 0;
  TheHMC.Parameters.MD.MDsteps = 100;
  TheHMC.Parameters.MD.trajL   = 1.0;
  
  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  
  TheHMC.Run();  // no smearing
  
  // Report final values
  RealD mean_action = solver.status().last_update.mean_action;
  RealD a_final  = constrained_action.parameters().a;
  std::cout << GridLogMessage << "Example_LLR: a_final = " << a_final << ", Sunconstrained = " << mean_action << std::endl;

  Grid_finalize();

} // main
