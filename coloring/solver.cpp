#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>


class Greedy {
public:
    Greedy(const std::string& fileName);
    ~Greedy() = default;

    Greedy(const Greedy&) = delete;
    Greedy(Greedy&&) = delete;
    Greedy& operator=(const Greedy&) = delete;
    Greedy& operator=(Greedy&&) = delete;

public:
    void run();
    std::pair<int, std::vector<int>> getAnswer() {
        return std::make_pair(colorIndex_, colors_);
    }

private:
    int getColor(int node);

private:
    static int UNCOLORED;

private:
    int numNode_;
    int numEdge_;
    int colorIndex_;

    std::vector<std::vector<int>> adjacency_;

    std::vector<int> colors_;
};


int Greedy::UNCOLORED = -1;


Greedy::Greedy(const std::string& fileName)
    :   colorIndex_(0) {

    std::ifstream input(fileName);
    input >> numNode_ >> numEdge_;

    adjacency_.resize(numNode_);
    colors_.resize(numNode_);

    for (int i = 0 ; i < numEdge_ ; ++i) {
        int src, dst;
        input >> src >> dst;
        adjacency_[src].push_back(dst);
        adjacency_[dst].push_back(src);
    }
}

void Greedy::run() {

    auto compare = [] (const std::pair<int, int>& lhs,
                       const std::pair<int, int>& rhs) {
        return lhs.second < rhs.second;
    };
    std::priority_queue<std::pair<int, int>,
                        std::vector<std::pair<int, int>>,
                        decltype(compare)> queue(compare);

    for (int i = 0 ; i < numNode_ ; ++i) {
        queue.push(std::make_pair(i, adjacency_[i].size()));
        colors_[i] = UNCOLORED;
    }

    while (!queue.empty()) {
        const auto& item = queue.top();
        int node = item.first;
        colors_[node] = getColor(node);
        queue.pop();
    }
}

int Greedy::getColor(int node) {

    if (colorIndex_ == 0) {
        return colorIndex_++;
    }

    std::vector<bool> occupation(colorIndex_);

    const auto& neighbors = adjacency_[node];
    for (auto neighbor : neighbors) {
        int color = colors_[neighbor];
        if (color == UNCOLORED) {
            continue;
        }

        occupation[color] = true;
    }

    for (int i = 0 ; i < colorIndex_ ; ++i) {
        if (occupation[i] == false) {
            return i;
        }
    }

    return colorIndex_++;
}


int main(int argc, char** argv) {

    if (argc != 2) {
        std::cerr << "Please specify the input data path" << std::endl;
        return -1;
    }

    std::string fileName(argv[1]);
    Greedy solver(fileName);
    solver.run();

    auto answer = solver.getAnswer();
    std::cout << answer.first << " 0" << std::endl;
    for (auto color : answer.second) {
        std::cout << color << ' ';
    }
    std::cout << std::endl;

    return 0;
}