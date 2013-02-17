/*!
 *	\brief			A PDE solver for the laplace's equation			
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

// Function prototypes
template <typename T>
void cmdLineArgToValue(const char * str, T & value, T defValue = T());

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
		
		value_type * begin(){ return &(m_data[0]); }
		value_type * end(){ return (m_size > 0)? &(m_data[m_size * m_size]) : 0; }
		
		std::size_t size(){ return m_size; }
	private:
		std::size_t m_size;	
		value_type * m_data;
};
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
	
	return 0;
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
