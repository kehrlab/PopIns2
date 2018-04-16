/*!
* \file    src/CompactedDBG_Data_extension.h
* \brief   Library for extending the functionality of the compacted DBG by implementing the data wrapper CDBG_Data_t
*
*/
#ifndef COMPACTED_DBG_DATA_EXTENSION_
#define COMPACTED_DBG_DATA_EXTENSION_



// =========================
// Includes
// =========================
#include <bifrost/CompactedDBG.hpp>

#include <cstdint>                      /* uint8_t */


// =========================
// Includes
// =========================
class DataExtension : public CDBG_Data_t<DataExtension> {

    public:

        DataExtension() : b(NOT_VISITED_SEEN), id(0) {}

        static void join(const UnitigMap<DataExtension>& um_dest, const UnitigMap<DataExtension>& um_src);
        static void sub(DataExtension* data_dest, const UnitigMap<DataExtension>& um_src, bool last_extraction);
        string serialize() const;

        unsigned getID() const {return id;}
        void setID(const unsigned id) {this->id = id;}

        inline void set_visited() { b = VISITED; } // Set the boolean to "visited"
        inline void set_seen() { b = SEEN; } // Set the boolean to "seen"
        inline void set_not_seen_visited() { b = NOT_VISITED_SEEN; } // Set the boolean to "not seen and not visited"

        inline bool is_visited() const { return (b == VISITED); } // return if the boolean is "visited"
        inline bool is_not_visited() const { return (b != VISITED); } // return if the boolean is "not visited"

        inline bool is_seen() const { return (b == SEEN); } // return if the boolean is "seen"
        inline bool is_not_seen() const { return (b != SEEN); } // return if the boolean is "not seen"

    private:

        const static uint8_t NOT_VISITED_SEEN = 0x0;
        const static uint8_t VISITED = 0x1;
        const static uint8_t SEEN = 0x2;

        uint8_t b;  // value of traversal status

        unsigned int id;
};





#endif /*COMPACTED_DBG_DATA_EXTENSION_*/