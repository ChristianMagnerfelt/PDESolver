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

int g_fpPrecision = 4;	//!< Floating point precision for streams
float g_e; // Maximum error ( epsilon )

/*!
 *	\brief	A matrix class implementing the RAII pinciple
 */
template <typename T>
class Matrix{
	public:
		typedef T value_type;
		
		Matrix(const Matrix &) = delete;				//!< Disable copying
		Matrix & operator=(const Matrix &) = delete;	//!< Disable copying
		
		Matrix (std::size_t size) : m_size(size), m_data(new value_type[m_size * m_size]){}
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
Matrix<T> & jacobi(Matrix<T> & gridG, Matrix<T> & gridT, std::size_t numIters, T & e);
template <typename T>
void printGrid(const Matrix<T> & grid);
template <typename T>
void initializeGrid(Matrix<T> & grid);
template <typename T>
T calculateDifference(Matrix<T> & gridG, Matrix<T> & gridT);


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
	
	// Initialize grid
	Matrix<float> gridG(g_gridSize + 2);
	Matrix<float> gridT(g_gridSize + 2);
	g_e = 0;	

	initializeGrid(gridG);
	initializeGrid(gridT);

	// Set precision for cout
	out << std::setprecision(g_fpPrecision) << std::fixed;

	double startTime, endTime;
	
	// Do calculations	
	startTime = omp_get_wtime();
	auto & result = jacobi(gridG, gridT, g_numIters, g_e);
	endTime = omp_get_wtime();
	std::cout << "Maximum error" << g_e << std::endl;
	std::cout << "Calculations took " << endTime - startTime << " seconds" << std::endl;
	
	// Print result
	printGrid(result);
	
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
Matrix<T> & jacobi(Matrix<T> & gridG, Matrix<T> & gridT, std::size_t numIters, T & e)
{
	for(std::size_t t = 0; t < numIters; ++t)
	{
		for(std::size_t i = 1; i < gridG.size() - 1; ++i)
		{
			for(std::size_t j = 1; j < gridG.size() - 1; ++j)
			{
						gridT[i][j] = (gridG[i][j-1] + gridG[i-1][j] + gridG[i+1][j] + gridG[i][j+1]) / static_cast<T>(4);
			}
		}
		e = calculateDifference(gridG, gridT);
		//std::cout << e << std::endl;
		Matrix<T>::swap(gridG, gridT);
	}
	return gridG;
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
 *	‪\brief	Calculates the maximum difference (epsilon) between two matrices
 */
template <typename T>
T calculateDifference(Matrix<T> & gridG, Matrix<T> & gridT)
{
	if(gridG.size() != gridT.size())
	{
		std::cerr << "Invalid matrix dimension" << std::endl;
		return -1.0f;
	}
	T e = 0;
	for(std::size_t i = 0; i < gridG.size(); ++i)
	{
		for(std::size_t j = 0; j < gridT.size(); ++j)
		{
			e = std::max(e, std::abs(gridG[i][j] - gridT[i][j]));		
		}
	}
	return e;
}
