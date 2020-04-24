/**
 * @file    src/debug_macros.h
 * @brief   Library for debug marcos in PopIns2
 */
#ifndef POPINS2_DEBUG_MACROS_
#define POPINS2_DEBUG_MACROS_


#ifdef DEBUG
#define DEBUG_PRINT_STATUS(msg) std::cout << msg << std::endl
#else
#define DEBUG_PRINT_STATUS(msg)
#endif


#ifdef DEBUG
#define DEBUG_PRINT_UCM_STATUS(msg) std::cout << get_unitig_id(ucm) << ": " << msg << std::endl
#else
#define DEBUG_PRINT_UCM_STATUS(msg)
#endif


#ifdef DEBUG
#define DEBUG_PRINT_PRE_STATUS(msg) std::cout << get_unitig_id(pre) << ": " << msg << std::endl
#else
#define DEBUG_PRINT_PRE_STATUS(msg)
#endif


#ifdef DEBUG
#define DEBUG_PRINT_SUC_STATUS(msg) std::cout << get_unitig_id(suc) << ": " << msg << std::endl
#else
#define DEBUG_PRINT_SUC_STATUS(msg)
#endif


#ifdef DEBUG
#define DEBUG_PRINT_PRE_LECC_BORDER std::cout << "Border: ID " << ue_pre->getID() << " - Kmer[" << pre.getMappedTail().rep().toString() << "]" << std::endl
#else
#define DEBUG_PRINT_PRE_LECC_BORDER
#endif


#ifdef DEBUG
#define DEBUG_PRINT_SUC_LECC_BORDER std::cout << "Border: ID " << ue_suc->getID() << " - Kmer[" << suc.getMappedHead().rep().toString() << "]" << std::endl
#else
#define DEBUG_PRINT_SUC_LECC_BORDER
#endif


#ifdef DEBUG
#define DEBUG_PRINT_PARTNER_TO_CHECK std::cout << "Border: ID " << data->getID() << " - Kmer[" << border2check.toString() << "]" << std::endl
#else
#define DEBUG_PRINT_PARTNER_TO_CHECK
#endif


#ifdef DEBUG
#define DEBUG_PRINT_TRACEBACK tb.print();
#else
#define DEBUG_PRINT_TRACEBACK
#endif



#endif /*POPINS2_DEBUG_MACROS_*/
