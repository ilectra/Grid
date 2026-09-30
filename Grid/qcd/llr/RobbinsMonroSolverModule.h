/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: RobbinsMonroSolverModule.h

Copyright (C) 2026

Author: Ryan Hill <Ryan.Hill@ed.ac.uk>

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

#pragma once
#include <string>

#include <Grid/qcd/llr/RobbinsMonroSolver.h>
#include <Grid/qcd/observables/hmc_observable.h>
#include <Grid/qcd/modules/ObservableModules.h>

NAMESPACE_BEGIN(Grid);

inline void LogRobbinsMonroUpdate(const RobbinsMonroUpdate &update,
                                  const std::string &context)
{
  std::cout << GridLogMessage << "[RobbinsMonro" << context << "] trajectory: " << update.trajectory
            << ", iteration: " << update.iteration
            << ", mean action: " << update.mean_action
            << ", residual: " << update.residual
            << ", a: " << update.previous_a << " -> " << update.updated_a << std::endl;
}

template <class RMSolver>
class LLRActionLogger : public HmcObservable<typename RMSolver::Field> {

public:
  typedef typename RMSolver::Field Field;
  
  LLRActionLogger(RMSolver &solver) : solver_(solver) {}

  virtual void TrajectoryComplete(int trajectory, Field &U,
                                  GridSerialRNG &, GridParallelRNG &)
  {
    int previous_iteration = solver_.status().iteration;
    solver_.record_configuration(trajectory, U);
    report_update(previous_iteration);
  }

  virtual void TrajectoryComplete(int trajectory,
                                  ConfigurationBase<Field> &configuration,
                                  GridSerialRNG &, GridParallelRNG &,
                                  bool)
  {
    int previous_iteration = solver_.status().iteration;
    solver_.record_configuration(trajectory, configuration);
    report_update(previous_iteration);
  }

  std::vector<RealD> getMeanActionVector()
  {
    //std::cout << "In getMeanActionVector " << mean_action_.size() << std::endl;
    std::vector<RealD> out;
    out = mean_action_;
    //std::cout << "In getMeanActionVector " << out.size() << std::endl;
    return out;
  }

private:
  RMSolver &solver_;
  std::vector<RealD> mean_action_;

  void report_update(int previous_iteration)
  {
    RobbinsMonroStatus status = solver_.status();
    if (status.iteration == previous_iteration)
    {
      return;
    }

    LogRobbinsMonroUpdate(status.last_update, "");
    //std::cout << mean_action_.size() << std::endl;
    mean_action_.push_back(status.last_update.mean_action);
    //std::cout << mean_action_.size() << std::endl;
    std::cout << GridLogMessage << "[Action Parameters]" << std::endl;
    std::cout << solver_.action_.LogParameters();
  }
};



template <class RMSolver>
class RobbinsMonroSolverModule
  : public ObservableModule<LLRActionLogger<RMSolver>, NoParameters> {
public:
  typedef typename RMSolver::Field Field;

  typedef ObservableModule<LLRActionLogger<RMSolver>, NoParameters> ObsBase;
  using ObsBase::ObsBase;

  RobbinsMonroSolverModule(RMSolver &solver): ObsBase(NoParameters()) {
    this->ObservablePtr.reset(new LLRActionLogger<RMSolver>(solver));
  }
  
  virtual void initialize(){}
};

NAMESPACE_END(Grid);

