#ifndef BIN_COUNT_H
#define BIN_COUNT_H

#include <vector>

namespace BinCount {
	/*
	* @brief This is just a binary counter.
	*/
	class Bin_count {
	public:
		/*
		* @brief Constructor
		*
		* @param [in] size The size of the counter.
		*/
		Bin_count(int size);

		/*
		* @brief Returns the current value of the counter at a given position.
		*
		* @param [in] pos Position to be returned.
		* @returns The current value of the counter at a given position.
		*/
		bool at(int pos);

		/*
		* @brief Increases the counter by one.
		*/
		void increase();

		/*
		* @brief Resets the counter to 0 in each bit.
		*/
		void reset();

		/*
		* @brief Prints the counter.
		*/
		void print();

		/*
		* @brief Checks if the counter can be increased.
		*
		* @returns Returns true if the counter cannot be increased anymore.
		*/
		bool is_full();

		/*
		* @brief Returns the size of the counter.
		*
		* @returns Size of the counter.
		*/
		int get_size();

	private:
		std::vector<bool> counter;
		int size;
		bool full;
	};
}
#endif
