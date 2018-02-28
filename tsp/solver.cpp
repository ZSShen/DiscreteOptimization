#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <unordered_set>
#include <cmath>
#include <time.h>


class TabuList {
public:
    TabuList(int capacity)
        : size_(0), capacity_(capacity) { }
    ~TabuList() = default;

    TabuList(const TabuList&) = delete;
    TabuList(TabuList&&) = delete;
    TabuList& operator=(const TabuList&) = delete;
    TabuList& operator=(TabuList&&) = delete;

public:
    bool find(int key) {
        auto iter = set_.find(key);
        return iter == set_.cend() ? false : true;
    }

    void insert(int key) {
        if (size_ == capacity_) {
            while (!queue_.empty()) {
                queue_.pop();
            }
            set_.clear();
            size_ = 0;
        }

        queue_.push(key);
        set_.insert(key);
        ++size_;
    }

private:
    int size_;
    int capacity_;

    std::priority_queue<int> queue_;
    std::unordered_multiset<int> set_;
};


class TwoOpt {
public:
    TwoOpt(const std::string& fileName,
           double errorRate,
           int dropRate,
           int maxRound);

    ~TwoOpt() = default;

    TwoOpt(const TwoOpt&) = delete;
    TwoOpt(TwoOpt&&) = delete;
    TwoOpt& operator=(const TwoOpt&) = delete;
    TwoOpt& operator=(TwoOpt&&) = delete;

public:
    void run();
    std::pair<double, std::vector<int>> getAnswer() {
        return make_pair(ansTourLen_, ansTour_);
    }

private:
    void swapEdge(int k, int h);

private:
    int numNode_;
    int numTrial_;

    double errorRate_;
    int dropRate_;
    int maxRound_;

    double ansTourLen_;
    double optTourLen_;
    double curTourLen_;

    std::vector<std::pair<double, double>> map_;
    std::vector<std::vector<double>> distance_;

    std::vector<int> ansTour_;
    std::vector<int> optTour_;
    std::vector<int> curTour_;
};


TwoOpt::TwoOpt(const std::string& fileName,
               double errorRate,
               int dropRate,
               int maxRound)
    :   numTrial_(1),
        errorRate_(errorRate),
        dropRate_(dropRate),
        maxRound_(maxRound) {

    std::ifstream input(fileName);
    input >> numNode_;

    map_.resize(numNode_);
    distance_.resize(numNode_);
    for (int i = 0 ; i < numNode_ ; ++i) {
        distance_[i].resize(numNode_);
        input >> map_[i].first >> map_[i].second;
    }

    for (int i = 0 ; i < numNode_ ; ++i) {
        for (int j = i + 1 ; j < numNode_ ; ++j) {
            double xSqr = pow(map_[j].first - map_[i].first, 2);
            double ySqr = pow(map_[j].second - map_[i].second, 2);
            distance_[i][j] = distance_[j][i] = sqrt(xSqr + ySqr);
        }
    }

    // Construct the initial tour
    optTourLen_ = 0;
    curTour_.resize(numNode_ + 1);
    optTour_.resize(numNode_ + 1);
    optTour_[0] = 0;

    for (int i = 1 ; i < numNode_ ; ++i) {
        optTour_[i] = i;
        optTourLen_ += distance_[i][i - 1];
    }

    optTour_[numNode_] = 0;
    optTourLen_ += distance_[0][numNode_ - 1];

    ansTourLen_ = optTourLen_;
}

void TwoOpt::run() {

    srand(time(NULL));

    int numFail = 0;
    int numRound = 0;
    TabuList tabu(numNode_ * 100);

    while (true) {
        // Continuously swap edges until no improvement can be found.
        for (int i = 1 ; i < numNode_ - 3 ; ++i) {
            for (int j = i + 1 ; j < numNode_ ; ++j) {

                ++numRound;
                if (numRound == maxRound_) {
                    return;
                }

                swapEdge(i, j);

                // Use Tabu list to gradually reduce the upper bound of
                // tour length found in each round.
                if (tabu.find(curTourLen_)) {
                    continue;
                }

                if (curTourLen_ < optTourLen_) {
                    optTourLen_ = curTourLen_;
                    optTour_ = curTour_;

                    if (curTourLen_ < ansTourLen_) {
                        ansTourLen_ = curTourLen_;
                        ansTour_ = curTour_;
                    }

#ifdef DEBUG
                    std::cout.precision(16);
                    std::cout << curTourLen_ << std::endl;
#endif

                    numFail = 0;
                    continue;
                }

                // If the current tour is not better than the local optimal,
                // we adopt Simulated Annealing to reheat the local search
                // with carefully crafted parameters.
                double diff = curTourLen_ - optTourLen_;
                if (diff / optTourLen_ < errorRate_) {
                    auto pick = rand() % dropRate_;
                    if (pick == 1) {
                        optTourLen_ = curTourLen_;
                        optTour_ = curTour_;
                        numFail = 0;
                    }
                    continue;
                }

                tabu.insert(curTourLen_);

                ++numFail;
                if (numFail == numNode_ * 10) {
                    return;
                }
            }
        }
    }
}

void TwoOpt::swapEdge(int k, int h) {

    for (int i = 0 ; i < k ; ++i) {
        curTour_[i] = optTour_[i];
    }

    for (int i = k, j = h ; i <= h ; ++i, --j) {
        curTour_[i] = optTour_[j];
    }

    for (int i = h + 1 ; i <= numNode_ ; ++i) {
        curTour_[i] = optTour_[i];
    }

    curTourLen_ = 0;
    for (int i = 1 ; i <= numNode_ ; ++i) {
        int src = curTour_[i - 1];
        int dst = curTour_[i];
        curTourLen_ += distance_[src][dst];
    }
}


int main(int argc, char** argv) {

    if (argc != 5) {
        std::cerr << "Please specify: "
                     "INPUT_PATH, "
                     "ERROR_RATE, "
                     "RANDOM_DROP_RATE, "
                     "MAX_ROUND" << std::endl;
        return -1;
    }

    std::string fileName(argv[1]);
    double errorRate = atof(argv[2]);
    int dropRate = atoi(argv[3]);
    int maxRound = atoi(argv[4]);

    TwoOpt solver(fileName, errorRate, dropRate, maxRound);
    solver.run();

    auto pair = solver.getAnswer();
    std::cout.precision(16);
    std::cout << pair.first << " 0" << std::endl;
    for (int i = 0 ; i < pair.second.size() - 1; ++i) {
        std::cout << pair.second[i] << ' ';
    }
    std::cout << std::endl;

    return 0;
}