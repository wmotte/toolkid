/**
 * Simple histogram implementation.
 */
class Histogram
{
private:

	std::vector< PixelType > histogram;
	unsigned int nBins;
	PixelType xLow;
	PixelType xHigh;
	PixelType delBin;

public:

	/**
	 * Construct histogram.
	 */
	Histogram( const std::vector< PixelType >& data, PixelType low, PixelType high, unsigned int bins )
	{
		xLow = low;
		xHigh = high;
		nBins = bins;

		histogram.assign( bins, 0 );

		delBin = ( xHigh - xLow ) / static_cast< PixelType > ( nBins );

		for ( unsigned int i = 0; i < data.size(); i++ )
		{
			PixelType dat = data[i];
			if ( dat < xLow )
				std::cerr << "underflow detected" << std::endl;

			else if ( dat >= xHigh )
				std::cerr << "overflow detected" << std::endl;

			else
			{
				unsigned int bin = static_cast< unsigned int > ( ( ( dat - xLow ) / delBin ) );
				if ( ( bin >= 0 ) && ( bin < nBins ) )
				{
					histogram[bin]++;
				}
			}
		}
	}

	/**
	 * Return histogram.
	 */
	std::vector< PixelType > getHistogram()
	{
		return histogram;
	}
};
