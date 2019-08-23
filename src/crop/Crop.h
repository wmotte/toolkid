
#ifndef CROP_H_
#define CROP_H_


template < int dimension >
class Crop {
public:
	void run ( const std::string& input, const std::string& output, float minValue );
};

#include "Crop.cpp"

#endif /* CROP_H_ */
