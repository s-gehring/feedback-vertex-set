#include <iostream>
#include "bin_count.hpp"

namespace BinCount {
	Bin_count::Bin_count(int size) : counter(size,0) {
		Bin_count::size = size;
		if (size > 0) {
			full = false;
		}
		else {
			std::cout << "Warning: Tried to create an empty counter!" << std::endl;
			full = true;
		}
	}

	bool Bin_count::at(int pos) {
		return counter[pos];
	}

	void Bin_count::increase() {
		if (!full) {
			int pos = size - 1;
			// first flip all ones
			while (pos >= 0 && counter[pos]) {
				counter[pos] = 0;
				pos--;
			}
			// then flip the rightmost zero (if the counter can be increased)
			if (pos < 0) {
				full = true;
			}
			else {
				counter[pos] = 1;
			}
		}
		else {
			std::cout << "Warning: Counter cannot be increased anymore. Please reset!" << std::endl;
		}
	} 

	void Bin_count::reset() {
		std::fill(counter.begin(), counter.end(), 0);
		full = false;
	}
    
    size_t Bin_count::get_true_values() {
        size_t x = 0;
        for(const auto &it:counter) {
            if(it) ++x;   
        }
        return x;
    }
    
    std::string Bin_count::to_string() const {
        std::string res = "";
        for(const auto &it : counter) {
            res = res + std::to_string(it) ;

        }
        return res;
    }
    
	void Bin_count::print() {
		std::cout<<to_string()<<std::endl;
	}

	bool Bin_count::is_full() {
		return full;
	}

	int Bin_count::get_size() {
		return size;
	}
}
