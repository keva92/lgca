/*
 * This file is part of LGCA, an implementation of a Lattice Gas Cellular Automaton
 * (https://github.com/keva92/lgca).
 *
 * Copyright (c) 2015-2017 Kerstin Vater, Niklas Kühl, Christian F. Janßen.
 *
 * LGCA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * LGCA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with lgca. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LGCA_BITSET_H_
#define LGCA_BITSET_H_

#include "lgca_common.h"

#include <tbb/spin_mutex.h>

#include <cstring> // memset

namespace lgca {

class Bitset
{
public:

    using Block = uint8_t;

    static constexpr Block BITS_PER_BLOCK = std::numeric_limits<Block>::digits;


    // A proxy class to simulate lvalues of bit type
    class reference
    {
        friend class Bitset;

        // The one and only non-copy constructor
        reference(Block& b, Block pos)
            : mBlock(b), mMask((assert(pos < BITS_PER_BLOCK), Block(1) << pos)) { }

        void operator&(); // Left undefined

    public:

        inline operator bool()  const { return (mBlock & mMask) != 0; }
        inline bool operator~() const { return (mBlock & mMask) == 0; }

        inline reference& flip() { /*tbb::spin_mutex::scoped_lock lock(mMutex);*/ do_flip(); return *this; }

        inline reference& operator=(bool x)               { do_assign(  x); return *this; } // For b[i] = x
        inline reference& operator=(const reference& rhs) { do_assign(rhs); return *this; } // For b[i] = b[j]

        inline reference& operator|=(bool x) { if  (x) do_set();   return *this; }
        inline reference& operator&=(bool x) { if (!x) do_reset(); return *this; }
        inline reference& operator^=(bool x) { if  (x) do_flip();  return *this; }
        inline reference& operator-=(bool x) { if  (x) do_reset(); return *this; }

     private:

        Block& mBlock;
        const Block mMask;

        inline void do_set() { mBlock |= mMask; }
        inline void do_reset() { mBlock &= ~mMask; }
        inline void do_flip() { mBlock ^= mMask; }
        inline void do_assign(bool x) { x? do_set() : do_reset(); }
    };

    Bitset() : m_bits(nullptr), m_num_bits(0), m_num_blocks(0) {}

    // Allocate a block of memory for an array of getNumBlocks(size) elements,
    // each of them sizeof(Block) bytes long, and initialize all its bits to zero
    Bitset(size_t size)
        : m_num_bits(size), m_num_blocks(get_num_blocks(size))
    {
        m_bits = (Block*)calloc(m_num_blocks, sizeof(Block));
    }

    ~Bitset()
    {
        free(m_bits);
        m_bits = nullptr;
    }

    // Resize bitset and set all its bits to zero
    inline void resize(size_t size)
    {
        if (m_bits) { free(m_bits); };

        m_num_bits   = size;
        m_num_blocks = get_num_blocks(size);

        m_bits = (Block*)calloc(m_num_blocks, sizeof(Block));
    }

    // Return the value of the bit at position pos
    inline bool operator[](size_t pos) const
    {
        assert(pos < m_num_bits);
        return (m_bits[block_idx(pos)] & bit_mask(pos)) != 0;
    }

    // Return a reference to the bit at position pos
    inline reference operator[](size_t pos)
    {
        assert(pos < m_num_bits);
        return reference(m_bits[block_idx(pos)], bit_idx(pos));
    }

    // Return the value of the block at position pos
    inline Block operator()(size_t pos) const
    {
        assert(pos < m_num_blocks);
        return m_bits[pos];
    }

    // Return a reference to the block at position pos
    inline Block& operator()(size_t pos)
    {
        assert(pos < m_num_blocks);
        return m_bits[pos];
    }

    // Set the bit at position pos to the value value
    inline void set(size_t pos, bool value = true)
    {
        assert(pos < m_num_bits);
        if (value) m_bits[block_idx(pos)] |= bit_mask(pos);
        else       reset(pos);
    }

    // Set all bits to true
    inline void set()
    {
        memset(m_bits, int(~Block(0)), m_num_blocks * sizeof(Block));
    }

    // Reset (to zero) the bit at position pos
    inline void reset(size_t pos)
    {
        assert(pos < m_num_bits);
        m_bits[block_idx(pos)] &= ~bit_mask(pos);
    }

    // Reset (to zero) all bits in the bitset
    inline void reset()
    {
        memset(m_bits, int(Block(0)), m_num_blocks * sizeof(Block));
    }

    // Flip the bit value at position pos (converting zeros into ones and ones into zeros)
    inline void flip(size_t pos)
    {
        assert(pos < m_num_bits);
        m_bits[block_idx(pos)] ^= bit_mask(pos);
    }

    // Flip all bit values in the bitset (converting zeros into ones and ones into zeros)
    inline void flip()
    {
        for (size_t i = 0; i < m_num_blocks; ++i)
            m_bits[i] = ~m_bits[i];
    }

    // TODO Return the number of bits in the bitset that are set (i.e., that have a value of one)
    inline size_t count() const
    {
        size_t counter = 0;
        for (size_t i = 0; i < m_num_bits; ++i)
            counter += this->operator [](i);

        return counter;
    }

    // Return the number of bits in the bitset
    inline size_t size() const
    {
        return m_num_bits;
    }

    inline void print() const
    {
        for (size_t i = 0; i < m_num_bits; ++i)
            cout << this->operator[](i) << " ";
        cout << endl;
    }

    inline void copy(const Bitset& other)
    {
        assert(other.size() == m_num_bits);
        memcpy(/*dst=*/(void*)m_bits, /*src=*/(const void*)other.ptr(), /*numBytes=*/m_num_blocks * sizeof(Block));
    }

    inline Block* ptr() { return &m_bits[0]; }
    inline const Block* ptr() const { return &m_bits[0]; }

    inline void operator=(Block* ptr)
    {
        // That's pretty dangerous... Let's hope bitset sizes do match...
        m_bits = ptr;
    }

    inline void fill_random()
    {
        for (size_t i = 0; i < m_num_bits; ++i)
            this->operator [](i) = (rand() % 2);
    }


private:

    static inline size_t block_idx(size_t pos) { return pos / BITS_PER_BLOCK; }
    static inline Block  bit_idx  (size_t pos) { return static_cast<Block>(pos % BITS_PER_BLOCK); }
    static inline Block  bit_mask (size_t pos) { return Block(1) << bit_idx(pos); }

    inline size_t get_num_blocks(size_t size)
    {
        return size / BITS_PER_BLOCK + static_cast<int>(size % BITS_PER_BLOCK != 0);
    }

    // Bits are represented as a linear array of Blocks, and the size of a Block is 8 bits
    Block* m_bits;

    size_t m_num_bits;
    size_t m_num_blocks;

}; // class Bitset

} // namespace lgca

#endif /* LGCA_BITSET_H_ */
