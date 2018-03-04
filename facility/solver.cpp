#include <ortools/base/commandlineflags.h>
#include <ortools/base/logging.h>
#include <ortools/linear_solver/linear_solver.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>


using namespace operations_research;


struct Facility {
    Facility(double setupCost, int capacity, double x, double y) {
        this->setupCost = setupCost;
        this->capacity = capacity;
        this->x = x;
        this->y = y;
    }

    double setupCost;
    int capacity;
    double x, y;
};


struct Client {
    Client(int demand, double x, double y) {
        this->demand = demand;
        this->x = x;
        this->y = y;
    }

    int demand;
    double x, y;
};


class IntegerProgramming {
public:
    IntegerProgramming(const std::string& fileName);
    ~IntegerProgramming() = default;

    IntegerProgramming(const IntegerProgramming&) = delete;
    IntegerProgramming(IntegerProgramming&&) = delete;

    IntegerProgramming& operator=(const IntegerProgramming&) = delete;
    IntegerProgramming& operator=(IntegerProgramming&&) = delete;

public:
    void run();

private:
    void run(MPSolver::OptimizationProblemType solverType);

private:
    int numFacility_;
    int numClient_;

    std::vector<Facility> facilities_;
    std::vector<Client> clients_;
};


IntegerProgramming::IntegerProgramming(const std::string& fileName) {

    std::ifstream input(fileName);
    input >> numFacility_ >> numClient_;

    for (int i = 0 ; i < numFacility_ ; ++i) {
        double setupCost;
        int capacity;
        double x, y;

        input >> setupCost >> capacity >> x >> y;
        facilities_.push_back(Facility(setupCost, capacity, x, y));
    }

    for (int i = 0 ; i < numClient_ ; ++i) {
        int demand;
        double x, y;

        input >> demand >> x >> y;
        clients_.push_back(Client(demand, x, y));
    }
}

void IntegerProgramming::run() {
    #if defined(USE_GLPK)
        run(MPSolver::GLPK_MIXED_INTEGER_PROGRAMMING);
    #endif
    #if defined(USE_CBC)
        run(MPSolver::CBC_MIXED_INTEGER_PROGRAMMING);
    #endif
    #if defined(USE_SCIP)
        run(MPSolver::SCIP_MIXED_INTEGER_PROGRAMMING);
    #endif
    #if defined(USE_GUROBI)
        run(MPSolver::GUROBI_MIXED_INTEGER_PROGRAMMING);
    #endif
    #if defined(USE_CPLEX)
        run(MPSolver::CPLEX_MIXED_INTEGER_PROGRAMMING);
    #endif
}

void IntegerProgramming::run(MPSolver::OptimizationProblemType solverType) {

    MPSolver solver("Facility", solverType);

    // Prepare MIP constants.
    std::vector<std::vector<double>> transportCToF(numClient_);
    for (int i = 0 ; i < numClient_ ; ++i) {
        for (int j = 0 ; j < numFacility_ ; ++j) {
            double sqrX = pow(clients_[i].x - facilities_[j].x, 2);
            double sqrY = pow(clients_[i].y - facilities_[j].y, 2);
            double dist = sqrt(sqrX + sqrY);
            transportCToF[i].push_back(dist);
        }
    }

    // Prepare MIP variables.
    std::vector<std::vector<MPVariable*>> assignCToF(numClient_);
    for (int i = 0 ; i < numClient_ ; ++i) {
        for (int j = 0 ; j < numFacility_ ; ++j) {
            std::string tag = "assign_c" +
                std::to_string(i) + "_f" + std::to_string(j);
            assignCToF[i].push_back(solver.MakeIntVar(0, 1, tag));
        }
    }

    std::vector<MPVariable*> openF(numFacility_);
    for (int i = 0 ; i < numFacility_ ; ++i) {
        std::string tag = "open_f" + std::to_string(i);
        openF[i] = solver.MakeIntVar(0, 1, tag);
    }

    // Prepare MIP optimal objective function.
    MPObjective* optima = solver.MutableObjective();
    for (int i = 0 ; i < numFacility_ ; ++i) {
        optima->SetCoefficient(openF[i], facilities_[i].setupCost);
    }

    for (int i = 0 ; i < numClient_ ; ++i) {
        for (int j = 0  ; j < numFacility_ ; ++j) {
            optima->SetCoefficient(assignCToF[i][j], transportCToF[i][j]);
        }
    }

    // Prepare MIP constraints.
    double infinity = solver.infinity();

    std::vector<std::vector<MPConstraint*>> openRelation(numClient_);
    for (int i = 0 ; i < numClient_ ; ++i) {
        for (int j = 0 ; i < numFacility_ ; ++j) {
            auto constraint = solver.MakeRowConstraint(-infinity, 0);
            constraint->SetCoefficient(assignCToF[i][j], 1);
            constraint->SetCoefficient(openF[j], -1);
            openRelation[i].push_back(constraint);
        }
    }

    std::vector<MPConstraint*> assignment(numClient_);
    for (int i = 0 ; i < numClient_ ; ++i) {
        auto constraint = solver.MakeRowConstraint(1, 1);
        for (int j = 0 ; i < numFacility_ ; ++j) {
            constraint->SetCoefficient(assignCToF[i][j], 1);
        }
        assignment[i] = constraint;
    }

    std::vector<MPConstraint*> capacity(numFacility_);
    for (int i = 0 ; i < numFacility_ ; ++i) {
        auto constraint = solver.MakeRowConstraint(
            -infinity, facilities_[i].capacity);
        for (int j = 0 ; j < numClient_ ; ++j) {
            constraint->SetCoefficient(assignCToF[j][i], clients_[j].demand);
        }
        capacity[i] = constraint;
    }

    //auto result = solver.Solve();
    //std::cout << optima->Value() << std::endl;
}


int main(int argc, char** argv) {

    if (argc != 2) {
        std::cerr << "Please specify the input data path" << std::endl;
        return -1;
    }

    std::string fileName(argv[1]);
    IntegerProgramming solver(fileName);

    solver.run();

    return 0;
}