//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BitChar - Fixed-size bitset of 8 bits
//

#ifndef BITCHAR_H
#define BITCHAR_H

struct BitChar
{
	public:
		BitChar() : d(0) {}
		
		// set the bit at idx to the value u
		void set(unsigned char idx, bool v);

		// returns true if bit at idx is set
		bool test(unsigned char idx) const;

		// flips the bit at idx
		void flip(unsigned idx);

		friend std::ostream& operator<<(std::ostream& out, const BitChar& bc);

	private:
		unsigned char d;
};

#endif
