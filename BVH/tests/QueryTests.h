#ifndef QUERY_TESTS_H
#define QUERY_TESTS_H


std::vector<double> TestPointToEdge(
        const std::vector<double>& p,
        const std::vector<std::vector<double>>& edgePts);

std::vector<double> TestPointToTri(
        const std::vector<double>& p,
        const std::vector<std::vector<double>>& triPts);


#endif
