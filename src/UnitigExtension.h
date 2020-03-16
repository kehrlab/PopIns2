/*!
* \file    src/ColoredCDBG_Data_extension.h
* \brief   Library for extending the functionality of the colored compacted DBG by implementing the data wrapper CCDBG_Data_t
*
*/
#ifndef UNITIG_EXTENSION_
#define UNITIG_EXTENSION_

// =========================
// Includes + defines
// =========================
#include <bifrost/ColoredCDBG.hpp>      /* has the CCDBG_Data_t template */
#include <bifrost/DataManager.hpp>
#include <bifrost/UnitigMap.hpp>

#include <string>                       /* std::to_string() */
#include <tuple>                        /* std::Pair */
#include <cstdint>                      /* uint8_t */


// =========================
// Structs
// =========================
/*!
* \class        UnitigExtension
* \headerfile   src/ColoredCDBG_Data_extension.h
* \brief        Struct extending the functionality of the colored compacted DBG.
* \details      The struct has to inherit from the CCDBG_Data_t struct and implement at least all its static functions.
* \ref          https://github.com/pmelsted/bfgraph/blob/master/src/ColoredCDBG.hpp
*/
struct UnitigExtension : public CCDBG_Data_t<UnitigExtension> {

    private:
        unsigned ID;
        float entropy;
        unsigned LECC;  // low entropy connected component

        const static uint8_t UNDISCOVERED = 0x0;
        const static uint8_t SEEN = 0x1;
        const static uint8_t VISITED = 0x2;

        uint8_t DFS_STATUS_FW;
        uint8_t DFS_STATUS_BW;

    public:

        // --------------
        // | Functions  |
        // --------------
        UnitigExtension() : ID(0),
                            entropy(-1),
                            LECC(0),
                            DFS_STATUS_FW(UNDISCOVERED),
                            DFS_STATUS_BW(UNDISCOVERED) {}

        unsigned getID() const {return ID;}
        void setID(const unsigned id) {ID = id;}

        float getEntropy() const {return entropy;}
        void setEntropy(const float e) {entropy = e;}

        float getLECC() const {return LECC;}
        void setLECC(const unsigned lecc) {LECC = lecc;}

        inline void set_undiscovered_fw() { DFS_STATUS_FW = UNDISCOVERED; }
        inline void set_seen_fw() { DFS_STATUS_FW = SEEN; }
        inline void set_visited_fw() { DFS_STATUS_FW = VISITED; }

        inline bool is_undiscovered_fw() const { return (DFS_STATUS_FW == UNDISCOVERED); }
        inline bool is_seen_fw() const { return (DFS_STATUS_FW == SEEN); }
        inline bool is_visited_fw() const { return (DFS_STATUS_FW == VISITED); }

        inline void set_undiscovered_bw() { DFS_STATUS_BW = UNDISCOVERED; }
        inline void set_seen_bw() { DFS_STATUS_BW = SEEN; }
        inline void set_visited_bw() { DFS_STATUS_BW = VISITED; }

        inline bool is_undiscovered_bw() const { return (DFS_STATUS_BW == UNDISCOVERED); }
        inline bool is_seen_bw() const { return (DFS_STATUS_BW == SEEN); }
        inline bool is_visited_bw() const { return (DFS_STATUS_BW == VISITED); }

        // -----------------------------------
        // | Implemented Abstract Functions  |
        // -----------------------------------
        /*
        void clear(const UnitigColorMap<UnitigExtension> &um_dest);

        void concat(const UnitigColorMap<UnitigExtension> &um_dest,
                    const UnitigColorMap<UnitigExtension> &um_src);

        void merge(const UnitigColorMap<UnitigExtension> &um_dest,
                   const const_UnitigColorMap<UnitigExtension> &um_src);

        void extract(const UnitigColors *uc_dest,
                     const UnitigColorMap<UnitigExtension> &um_src,
                     const bool last_extraction);
        */

        string serialize(const const_UnitigColorMap<UnitigExtension> &um_src) const;
};



#endif /*UNITIG_EXTENSION_*/
