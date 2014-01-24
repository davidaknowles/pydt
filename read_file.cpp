#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <iterator>
#include <boost/multi_array.hpp>

using namespace std;

std::vector<double> readRow(std::string row) {
    std::vector<double> retval;
    std::istringstream is(row);
    double num;
    while (is >> num)
        retval.push_back(num);

    return retval;
}

std::vector<std::vector<double> > readVector(std::istream &is) {
    std::string line;
    std::vector<std::vector<double> > retval;
    while (std::getline(is, line))
        retval.push_back(readRow(line));
    return retval;
}

std::vector<std::vector<double> > readFile(char *filename)
{
    ifstream in(filename);
    std::vector<std::vector<double> > retval = readVector(in);
    in.close();
    return retval;
}

std::vector<string > readNames(char *filename)
{
    std::vector<string > retval;
    ifstream in(filename);
    std::string line;

    std::getline(in, line);
    std::istringstream is(line);
    string st;
    while (is >> st)
        retval.push_back(st);
    in.close();
    return retval;
}


int main (int argc,char *argv[]) {

//   std::vector<std::vector<double> > r = readFile("sensitivity.txt");
//   for (int i=0;i<r.size();i++)
//   {
//     cout << r[i].size();
//     /*for (int j=0;j<r[i].size();j++)
//     {
//       //cout << r[i][j] << " ";
//
//     */
//     cout << endl;
//    }

    std::vector<string > r = readNames("cell_names.txt");
    for (int i=0; i<r.size(); i++)
        cout << r[i] << endl;

    return 0;
}
