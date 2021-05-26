#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

// function to read txt file and adding values to vector elements
Eigen::VectorXf CSVtoVector(std::string file, int N) {

    double data[N+1];

    Eigen::VectorXf vector = Eigen::VectorXf::Zero(N+1);

    std::ifstream input(file);

    for (int j = 0; j < N+1; ++j) {
        input >> data[j];
        vector[j] = data[j];

    }
    
    return vector;
}

// Saving the w vector to the file
void WriteFile(Eigen::VectorXf vector, std::string filename) {
    
    std::ofstream myfile (filename.append(".txt"), std::ios_base::app);

    if (myfile.is_open()){

            myfile << vector << "\n";
    }
    myfile.close();
}
    
// Performing LDU decomposition
std::vector<Eigen::MatrixXf> LUdecomposition(Eigen::MatrixXf M){

    std::vector<Eigen::MatrixXf> myVec;

    Eigen::MatrixXf U = Eigen::MatrixXf::Zero(M.rows(), M.cols());
    U.triangularView<Eigen::StrictlyUpper>() = M.triangularView<Eigen::StrictlyUpper>();
    
    Eigen::MatrixXf L = Eigen::MatrixXf::Zero(M.rows(), M.cols());
    L.triangularView<Eigen::StrictlyLower>() = M.triangularView<Eigen::StrictlyLower>();
    
    Eigen::MatrixXf D = Eigen::MatrixXf::Zero(M.rows(), M.cols());
    for(int i=0; i < M.rows(); ++i){

        D(i, i) = M(i, i);

    }
    
    myVec.push_back(L);
    myVec.push_back(D);
    myVec.push_back(U);

    return myVec;
}
// Performing LDU decomposition ends

// Building vector f
Eigen::VectorXf BoundaryConditions(Eigen::VectorXf f, Eigen::VectorXf Y, Eigen::VectorXf Z, Eigen::VectorXf P, int i, int n, float delt, float beta){

    f.block(0, 0, i, 1) = Y.block(0, 0, i, 1);
        
    for(int j=0; j < n; ++j){

        f.block(2*i*j + i, 0, i, 1) = -1.0 * delt * Z.block(i*j, 0, i, 1);

    }
    
    f.block(2*i*n + i, 0, i, 1) = -1.0 * beta * delt * Z.block(i*n, 0, i, 1);

    return f;
}
// Building vector f ends

// Building Matrix A
Eigen::MatrixXf BuildA(int i){

    float delx = 1.0/i;

    Eigen::MatrixXf A =  Eigen::MatrixXf::Identity(i+1, i+1);

    for(int j = 1; j < i; ++j){

        A(j, j) = 1.0/(delx*delx) * 2.0;

    }

    for(int j=1; j < i-1; ++j){

        A(j+1, j) = 1.0/(delx*delx) * -1.0;

    }

    for(int j=1; j < i-1; ++j){

        A(j, j+1) = 1.0/(delx*delx) * -1.0;

    }

    return A;
} 
// Building Matrix A ends

// Building Jacobi Smoother
Eigen::VectorXf JacobiSmoother(Eigen::MatrixXf A, Eigen::VectorXf w_old, Eigen::VectorXf f, int Iterations, float Tol){

    Eigen::VectorXf d;
    Eigen::VectorXf w_new;
    Eigen::VectorXf w_diff;

    float norm = 1.0;
    float omega = 1.0;

    std::vector<Eigen::MatrixXf> myVec= LUdecomposition(A);
    Eigen::MatrixXf L = myVec[0];
    Eigen::MatrixXf D = myVec[1];
    Eigen::MatrixXf U = myVec[2];

    // adding initial line to the file
    std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << "Residual: " << "\n";
        }
        myfile.close();

    // adding initial line to the file
    std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << "Iterations: " << "\n";
        }
        myfile_2.close();    

    int j = 0;
    while( (norm > Tol) && (j < Iterations) ){
        
        // calculating the defect
        d = f - A * w_old;

        d = D.inverse() * d;

        // Update the vector
        w_new = w_old + omega * d;

        // calculating the norm
        w_diff = w_new - w_old;

        //std::cout << "\n w_diff: \n" << w_diff;

        norm = w_diff.norm();

        // Printing the residual norm
        if(j==0){
            std::cout << "\n Residual norm: " << "\n" <<  std::endl;
        }
        std::cout << norm << "\n" <<  std::endl;
        
        // writing residual norm to the file
        std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << norm << "\n";
        }
        myfile.close();

        // writing number of iterations to the file
        std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << j << "\n";
        }
        myfile_2.close();

        w_old = w_new;
       
        ++j;

    }
    std::cout << "\n Solution Converged in: "  << j << " Iterations." << std::endl;

    return w_new;
}
// Building Jacobi Smoother ends

// Restriction operation starts here
Eigen::MatrixXf Restriction_operator(int i){
        
    Eigen::MatrixXf restricted_operator = Eigen::MatrixXf::Zero((i/2+1), i+1);
    
    Eigen::MatrixXf d_m(1, 3);

    d_m << 1, 2, 1;
    
    for(int j=1; j < i/2; ++j){

        restricted_operator.block(j, 2*j - 1, 1, 3) = (1.0/4) * d_m;
    }

    restricted_operator(0, 0) = 1;
    restricted_operator(i/2, i) = 1;

    return restricted_operator;
}
// Restriction operation ends here

// Restriction operation starts here
Eigen::VectorXf Restriction(Eigen::VectorXf d, Eigen::MatrixXf R_s){

    Eigen::VectorXf d_restricted = R_s * d;
        
    return d_restricted;
}
// Restriction operation ends here

// Restriction operation starts here
Eigen::VectorXf Prolongation(Eigen::VectorXf w, Eigen::MatrixXf P_s){

    Eigen::VectorXf w_prolongated = P_s * w;
        
    return w_prolongated;
}   
// Restriction operation ends here

// V-Cycle starts here
Eigen::VectorXf VCycle(Eigen::VectorXf w_old, Eigen::VectorXf f, int N){

    // V-Cycle starts
    Eigen::MatrixXf A = BuildA(N);

    Eigen::VectorXf w_tilda = JacobiSmoother(A, w_old, f, 20000, 1E-01);
    
    WriteFile(w_tilda, "w_tilda");

    Eigen::MatrixXf R_s = Restriction_operator(N);

    Eigen::MatrixXf P_s = 2 * R_s.transpose();
    P_s(0, 0) = 1;
    P_s(N, N/2) = 1;

    Eigen::VectorXf d = f - A * w_tilda;

    std::cout << "\n \n Restricting.... \n" << std::endl;

    Eigen::VectorXf d_tilda = Restriction(d, R_s); 

    Eigen::MatrixXf A_2h = R_s * A * P_s;

    // Coarse grid equation
    Eigen::VectorXf d_tilda_tilda = JacobiSmoother(A_2h, Eigen::VectorXf::Zero(d_tilda.size()), d_tilda, 20000, 1E-08);

    WriteFile(d_tilda_tilda, "d_tilda_tilda");

    std::cout << "\n \n Prolongating.... \n" << std::endl;

    Eigen::VectorXf d_tilda_tilda_vector = Prolongation(d_tilda_tilda, P_s);

    WriteFile(d_tilda_tilda_vector, "d_tilda_tilda_vector");

    Eigen::VectorXf w_tilda_vector = w_tilda + d_tilda_tilda_vector; 

    WriteFile(w_tilda_vector, "w_tilda_vector");
    
    Eigen::VectorXf w = JacobiSmoother(A, w_tilda_vector, f, 20000, 1E-04); 

    WriteFile(w, "w");
    
    return w;
}   
// V-cycle ends

// Single Cycle starts here
Eigen::VectorXf SingleCycle(Eigen::VectorXf w_old, Eigen::VectorXf f, int N){

    // V-Cycle starts
    Eigen::MatrixXf A = BuildA(N);

    Eigen::VectorXf w = JacobiSmoother(A, w_old, f, 50000, 1E-04);

    WriteFile(w, "w");

    return w;
}   
// Single cycle ends

int main(){

    remove("Residual.txt");
    remove("Iterations.txt");
    remove("W.txt");

    int N;
    
    std::cout << "Enter the number of space grids: " << std::endl;

    std::cin >> N;

    std::cout << "\n delx: " << 1.0/N << std::endl;

    Eigen::VectorXf Y = CSVtoVector("Y.txt", N);
    Eigen::VectorXf f = CSVtoVector("F.txt", N);

    std::cout << "\n Y: \n" << Y << std::endl;

    f(0) = Y(0);
    f(N) = Y(N);

    std::cout << "\n f: \n" << f << std::endl;

    Eigen::VectorXf w_old = Eigen::VectorXf::Random(N+1);

    //Eigen::VectorXf w1 = SingleCycle(w_old, f, N);

    Eigen::VectorXf w2 = VCycle(w_old, f, N);

    return 0;
}  