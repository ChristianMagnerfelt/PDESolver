/*!
 *	\brief			A PDE solver for the laplace's equationmake			
 *	\description	
 *	\author			Christian Magnerfelt
 *	\email			magnerf@kth.se
 */
#include <omp.h>

#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <algorithm>

#include <cmath>

// Global variables
const std::size_t DEFAULT_GRID_SIZE = 10;
const std::size_t DEFAULT_NUM_ITERS = 10;
const std::size_t DEFAULT_NUM_WORKERS = 10;
const std::string DEFAULT_OUT_FILENAME ("result.out");

std::size_t g_gridSize = DEFAULT_GRID_SIZE;
std::size_t g_numIters = DEFAULT_NUM_ITERS;
std::size_t g_numWorkers = DEFAULT_NUM_WORKERS;
std::string g_outFileName = DEFAULT_OUT_FILENAME;

std::fstream out;

int g_dfpPrecision = 15;	//!< Double floating point precision for streams

/*!
 *	\brief	A matrix class implementing the RAII pinciple
 */
template <typename T>
class Matrix{
	public:
		typedef T value_type;
		
		Matrix(const Matrix &) = delete;				//!< Disable copying
		Matrix & operator=(const Matrix &) = delete;	//!< Disable copying
		
		explicit Matrix (std::size_t size) 
			: m_size(size), m_data(new value_type[m_size * m_size]){}
		Matrix (std::size_t size, const T & defValue) 
			: m_size(size), m_data(new value_type[m_size * m_size])
		{std::fill(begin(), end(), defValue);}

		~Matrix(){ delete [] m_data; }
		
		value_type * operator [](std::size_t index){ return &(m_data[m_size * index]); }
		const value_type * operator [](std::size_t index) const { return &(m_data[m_size * index]); }
		
		value_type * begin(){ return &(m_data[0]); }
		value_type * end(){ return (m_size > 0)? &(m_data[m_size * m_size]) : 0; }
		
		std::size_t size() const { return m_size; }
		
		static void swap(Matrix<T> & matA, Matrix<T> & matB);
	private:
		std::size_t m_size;	
		value_type * m_data;
};

// Matrix non-member functions ( Grid specific )
template <typename T> 
Matrix<T> & jacobi(Matrix<T> & gridG, Matrix<T> & gridT);
template <typename T>
void printGrid(const Matrix<T> & grid);
template <typename T>
void initializeGrid(Matrix<T> & grid);
template <typename T>
T calculateMaxDifference(Matrix<T> & gridG, Matrix<T> & gridT);


// Function prototypes
template <typename T>
void cmdLineArgToValue(const char * str, T & value, T defValue = T());


/*!
 *	\brief	Program entry point
 */
int main(int argc, const char * argv [])
{
	// Read and extract command line arguments
	if(argc > 1)
	{
		switch(argc)
		{	
			default :
			case (5) :
				cmdLineArgToValue(argv[4], g_outFileName, DEFAULT_OUT_FILENAME);
			case (4) :
  				cmdLineArgToValue(argv[3], g_numWorkers, DEFAULT_NUM_WORKERS);
			case (3) :
				cmdLineArgToValue(argv[2], g_numIters, DEFAULT_NUM_ITERS);			
			case (2) :	
				cmdLineArgToValue(argv[1], g_gridSize, DEFAULT_GRID_SIZE);		
		}
	}
	out.open(g_outFileName.data(), std::ios_base::out);
	if(out.fail())
	{
		std::cerr << "Failed to open file " << g_outFileName.data() << std::endl;
		return 0;
	
	}
	std::cout << "Grid size is " << g_gridSize << std::endl;
	std::cout << "Number of iterations is " << g_numIters << std::endl;
	std::cout << "Number of workers is " << g_numWorkers << std::endl;
	
	// Allocate and initialize grid
	Matrix<double> gridG(g_gridSize + 2);
	Matrix<double> gridT(g_gridSize + 2);
	Matrix<double> compareGrid(g_gridSize + 2, 1.0);

	initializeGrid(gridG);
	initializeGrid(gridT);

	// Initialize other values
	double maxDiff = 0;			// Maximum difference between two iterations
	double maxError = 0;		// Maximum error between result and correct answer
	double startTime = 0;
	double endTime = 0;
	std::size_t totalIters = 0;
	std::size_t checkMaxDiffFactor = 2 * g_gridSize; // How often we check the maximum difference between two grids.

	// Set precision for cout
	out << std::setprecision(g_dfpPrecision) << std::fixed;

	// Set number of workers
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(g_numWorkers); // Use 4 threads for all consecutive parallel regions

	// Do calculations	
	startTime = omp_get_wtime();
	for(std::size_t t = 0; t < g_numIters; ++t)
	{
		jacobi(gridG, gridT);
		Matrix<double>::swap(gridG, gridT);

		// If the grid size is large its likely to converge slower, thus we do 
		// not need to check the maximum difference that often.
		if((t % checkMaxDiffFactor) == 0)
		{
			maxDiff = calculateMaxDifference(gridT, gridG);
			// If the maximum difference between the grid are 0 then break the loop 
			// as there are not further improvements to make. We are done!
			if(maxDiff <= 0.0)
			{
				totalIters = t + 1;
				break;
			}
		}
	}	
	maxError = calculateMaxDifference(gridG, compareGrid);
	endTime = omp_get_wtime();
	std::cout << "Completed " << totalIters << "/" << g_numIters << " iterations" << std::endl; 
	std::cout << "Maximum difference " << maxDiff << std::endl;
	std::cout << "Maximum error " << maxError << std::endl;
	std::cout << "Calculations took " << endTime - startTime << " seconds" << std::endl;
	
	// Print result
	printGrid(gridG);
	
	out.close();
	return 0;
}
/*!
 *	\brief	Initializes grid border values
 */
template <typename T>
void initializeGrid(Matrix<T> & grid)
{
	if(grid.size() < 3)
		return;
		
	// Fill first row
	for(std::size_t i = 0; i < grid.size(); ++i)
		grid[0][i] = static_cast<T>(1);

	// Fill last row
	for(std::size_t i = 0; i < grid.size(); ++i)
		grid[grid.size() - 1][i] = static_cast<T>(1);
		
	// Fill first column
	for(std::size_t i = 1; i < grid.size() - 1; ++i)
		grid[i][0] = static_cast<T>(1);
		
	// Fill last column
	for(std::size_t i = 1; i < grid.size() - 1; ++i)
		grid[i][grid.size() - 1] = static_cast<T>(1);
}
/*!
 *	\brief	Does jacobi iteration on grid
 */
template <typename T> 
Matrix<T> & jacobi(Matrix<T> & gridG, Matrix<T> & gridT)
{
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for(std::size_t i = 1; i < gridG.size() - 1; ++i)
		{
			for(std::size_t j = 1; j < gridG.size() - 1; ++j)
			{
				gridT[i][j] = (gridG[i][j-1] + gridG[i-1][j] + gridG[i+1][j] + gridG[i][j+1]) / static_cast<T>(4);
			}
		}
	}
	return gridT;
}
/*!
 *	\brief	Prints grid to cout
 */
template <typename T>
void printGrid(const Matrix<T> & grid)
{
	for(std::size_t i = 0; i < grid.size(); ++i)
	{
		for(std::size_t j = 0; j < grid.size(); ++j)
		{
			out << grid[i][j] << " ";
		}
		out << std::endl;
	}
}
/*!	
 *	\brief	Takes a command line argument and tries to stream to a value. If streaming fails
 * 			then value is assigned a default value.
 */
template <typename T> 
void cmdLineArgToValue(const char * str, T & value, T defValue = T())
{
	std::istringstream iss(str);
	iss >> value;
	if(iss.fail())
	{
		value = defValue;
		std::cerr << "String: " << str << " is not a valid format" << std::endl;
	}
}
/*!
 *	\brief	Swaps internal pointers of two matrices
 *
 */
template <typename T>
void Matrix<T>::swap(Matrix<T> & matA, Matrix<T> & matB)
{
	value_type * tmp = matA.m_data;
	matA.m_data = matB.m_data;
	matB.m_data = tmp;
}
/*!
 *	‪\brief	Calculates the maximum difference between two matrices
 */
template <typename T>
T calculateMaxDifference(Matrix<T> & gridG, Matrix<T> & gridT)
{
	if(gridG.size() != gridT.size())
	{
		std::cerr << "Invalid matrix dimension" << std::endl;
		return -1.0f;
	}
	T maxDiff = 0;
	for(std::size_t i = 0; i < gridG.size(); ++i)
	{
		for(std::size_t j = 0; j < gridT.size(); ++j)
		{
			maxDiff = std::max(maxDiff, std::abs(gridG[i][j] - gridT[i][j]));		
		}
	}
	return maxDiff;
}
