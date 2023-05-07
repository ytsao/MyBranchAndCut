#include <iostream>
#include "ilconcert/cplexconcertdoc.h"
#include "ilcplex/ilocplex.h"

ILOSTLBEGIN

// Branch on var with largest objective coefficient
// among those with largest infeasibility
ILOCPLEXGOAL1(MyBranchGoal, IloNumVarArray, vars) {
    IloNumArray x;
    IloNumArray obj;
    IntegerFeasibilityArray feas;

    x = IloNumArray(getEnv());
    obj = IloNumArray(getEnv());
    feas = IntegerFeasibilityArray(getEnv());
    getValues(x, vars);
    getObjCoefs(obj, vars);
    getFeasibilities(feas, vars);

    IloInt bestj = -1;
    IloNum maxinf = 0.0;
    IloNum maxobj = 0.0;
    IloInt cols = vars.getSize();
    for (IloInt j = 0; j < cols; j++) 
    {
        if (feas[j] == Infeasible) 
        {
            IloNum xj_inf = x[j] - IloFloor(x[j]);
            if (xj_inf > 0.5)
                xj_inf = 1.0 - xj_inf;
            if (xj_inf >= maxinf && 
                (xj_inf > maxinf || IloAbs(obj[j]) >= maxobj)) 
            {
                bestj = j;
                maxinf = xj_inf;
                maxobj = IloAbs(obj[j]);
            }
        }
    }

    IloCplex::Goal res;
    if (bestj >= 0) 
    {
        res = AndGoal(OrGoal(vars[bestj] >= IloFloor(x[bestj]) + 1,
            vars[bestj] <= IloFloor(x[bestj])),
            this);
    }

    x.end();
    obj.end();
    feas.end();

    return res;
}


int main()
{
    try 
    {
        IloEnv env;
        IloModel model(env);
        IloCplex cplex(model);

        IloObjective obj;
        IloNumVarArray var(env);
        IloRangeArray rng(env);

        // Import lp file into CPLEX
        cplex.importModel(model, "location.lp", obj, var, rng);
        cplex.extract(model);
        cplex.setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional); // Cplex would switch to Traditional B&C automaticly
#pragma warning(disable:4996)
        cplex.setParam(IloCplex::Param::Threads, atoi(getenv("NUMBER_OF_PROCESSORS")));

        cplex.solve(MyBranchGoal(env,var));

        env.out() << "Objective Value = " << cplex.getObjValue() << "\n";
        env.out() << "Finished\n";
    }
    catch (IloException ex) 
    {
        std::cout << ex.getMessage() << "\n";
    }
}
