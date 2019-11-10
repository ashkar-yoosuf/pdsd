#include <iostream>
#include <fstream>
#include <vector>
#include <sstream> //to use stringsteam
#include <algorithm>

using namespace std;
ifstream inFile;

int find_density(vector<int> answer)
{
    int degree = 0;
    ifstream input_file;
    input_file.open("./edges.txt");

    string node1, node2;

    int src = 0, dest = 0;

    int line = 0;

    while (input_file.good())
    {

        //for csv
        // getline(input_file, node1, ',');

        getline(input_file, node1, ' ');
        getline(input_file, node2, '\n');

        stringstream start_node(node1), end_node(node2);

        start_node >> src;
        end_node >> dest;

        if (find(answer.begin(), answer.end(), src) != answer.end() && find(answer.begin(), answer.end(), dest) != answer.end())
        {
            degree += 2;
        }
    }

    cout << (degree / (2.0 * (answer.size()))) << "\n";
}

vector<int> Find_Densest_Subgraph(int num_of_nodes, int num_of_edges)
{

    int min_degree = 0;
    int max_degree = num_of_edges;
    vector<int> subgraph;
    int difference = 1.0 / (num_of_nodes * (num_of_nodes - 1));

    while (max_degree - min_degree >= difference)
    {

        cout << "Algo at work!";
        int least_density = (max_degree + min_degree) / 2.0;

        vector<int> source_segment = make_graph(num_of_nodes, num_of_edges, least_density);

        if (!source_segment.empty()) //check if the vector is empty or not
        {
            min_degree = least_density;
            subgraph = source_segment;
        }
        else
        {
            max_degree = least_density;
        }

        return subgraph;
    }

    vector<int> make_graph(int num_of_nodes, int num_of_edges, int least_density)
    {

        return vector<int>{
            1,
            2,
        };
    }
}

int main()
{
    vector<int> numbers = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    find_density(numbers);
}
