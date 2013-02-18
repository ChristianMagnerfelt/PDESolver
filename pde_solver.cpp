/*!
 *	\brief			A PDE solver for the laplace's equationmake			
 *	\description	
 *	\author			Christian Magnerfelt
 *	\email			magnerf@kth.se
 */

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

// Global variables
const std::size_t DEFAULT_GRID_SIZE = 10;
const std::size_t DEFAULT_NUM_ITERS = 10;
const std::size_t DEFAULT_NUM_WORKERS = 10;

std::size_t g_gridSize = DEFAULT_GRID_SIZE;
std::size_t g_numIters = DEFAULT_NUM_ITERS;
std::size_t g_numWorkers = DEFAULT_NUM_WORKERS;

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
	private:
		std::size_t m_size;	
		value_type * m_data;
};

// Function prototypes
template <typename T>
void cmdLineArgToValue(const char * str, T & value, T defValue = T());
template <typename T>
void initializeGrid(Matrix<T> & grid);
template <typename T> 
Matrix<T> & jacobi(Matrix<T> & gridG, Matrix<T> & gridT, std::size_t numIters);
template <typename T>
void printGrid(const Matrix<T> & grid);

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
			case (4) :
  				cmdLineArgToValue(argv[3], g_numWorkers, DEFAULT_NUM_WORKERS);
			case (3) :
				cmdLineArgToValue(argv[2], g_numIters, DEFAULT_NUM_ITERS);			
			case (2) :	
				cmdLineArgToValue(argv[1], g_gridSize, DEFAULT_GRID_SIZE);		
		}
	}
	std::cout << "Grid size is " << g_gridSize << std::endl;
	std::cout << "Number of iterations is " << g_numIters << std::endl;
	std::cout << "Number of workers is " << g_numWorkers << std::endl;
	
	// Initialize grid
	Matrix<float> gridG(g_gridSize + 2);
	Matrix<float> gridT(g_gridSize + 2);
	
	initializeGrid(gridG);
	initializeGrid(gridT);
	
	auto & result = jacobi(gridG, gridT, g_numIters);
	
	printGrid(result);
	
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
Matrix<T> & jacobi(Matrix<T> & gridG, Matrix<T> & gridT, std::size_t numIters)
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
			std::cout << grid[i][j] << " ";
		}
		std::cout << std::endl;
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
