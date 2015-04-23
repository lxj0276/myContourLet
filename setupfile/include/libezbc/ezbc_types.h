#ifndef __EZBC_TYPES_H
#define __EZBC_TYPES_H

#include "stdint.h"

typedef uint8_t ezbc_byte_t;
typedef uint16_t ezbc_coord_t;
typedef uint32_t ezbc_point_t;
typedef float ezbc_real_t;

/* coefficient representation: */
/* sign   magnitude  padding   length  */
/*  31  30    30-length       4      0 */
/*   S   M ..... m  0 0 0 0 0 L .... l */
/* length = 0 => FLC */
/* length = 1-26 => VLC */
#define EZBC_LSB               0
#define EZBC_MSB              30
#define EZBC_SIGN      (1UL << 31)          /* sign of the codeword */
#define EZBC_LENGTH  (~(-1  << 5))          /* length of the codeword */
#define EZBC_ABS ((~EZBC_SIGN)-EZBC_LENGTH) /* magnitude of the codeword */

typedef uint32_t ezbc_coeff_t;

#define ezbc_coord_x(v) ((v) & 0xffff)
#define ezbc_coord_y(v) ((v) >> 16)
#define ezbc_point(x, y) (((y) << 16) | (x))

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

#define shift_ceil(x, l) (((x) + (1 << (l)) - 1) >> (l))

#define EZBC_CONTEXT_SIZE 2

#define TL_POS    0 /* Bit-pos for significance of top-left neighbour        */
#define TC_POS    1 /* Bit-pos for significance of top-centre neighbour      */
#define TC2_POS   2 /*2nd Bit_pos for significance of top-centre neighbour   */
#define TR_POS    3 /* Bit-pos for significance of top-right neighbour       */
#define CL_POS    4 /* Bit-pos for significance of centre-left neighbour     */
#define CL2_POS   5 /*2nd Bit-pos for significance of centre-left neighbour  */
#define CR_POS    6 /* Bit-pos for significance of centre-right neighbour    */
#define CR2_POS   7 /*2nd Bit-pos for significance of centre-right neighbour */
#define BL_POS    8 /* Bit-pos for significance of bottom-left neighbour     */
#define BC_POS    9 /* Bit-pos for significance of bottom-centre neighbour   */
#define BC2_POS  10 /*2nd Bit-pos for significance of bottom-centre neighbour*/
#define BR_POS   11 /* Bit-pos for significance of bottom-right neighbour    */
#define PA_POS   12 /* Bit-pos for significance of parent pel                */
                                                                             
#define OUT_OF_BOUNDS_POS 15 /* May be used to identify context words which
                                lie beyond the boundaries of the code block */

#define TL_SIG        (1<<TL_POS)
#define TC_SIG        (1<<TC_POS)
#define TR_SIG        (1<<TR_POS)
#define CL_SIG        (1<<CL_POS)
#define CR_SIG        (1<<CR_POS)
#define BL_SIG        (1<<BL_POS)
#define BC_SIG        (1<<BC_POS)
#define BR_SIG        (1<<BR_POS)
#define TC2_SIG       (1<<TC2_POS)
#define CL2_SIG       (1<<CL2_POS)
#define CR2_SIG       (1<<CR2_POS)
#define BC2_SIG       (1<<BC2_POS)
#define PA_SIG        (1<<PA_POS)
#define OUT_OF_BOUNDS (1<<OUT_OF_BOUNDS_POS)

/* zero coding context */
#define ZC_BITS 13
#define ZC_MASK (~((-1) << ZC_BITS))

#define V_PVE_BIT_POS 0
#define V_NVE_BIT_POS 1
#define H_PVE_BIT_POS 2
#define H_NVE_BIT_POS 3
#define NW_PVE_BIT_POS 4
#define NW_NVE_BIT_POS 5
#define NE_PVE_BIT_POS 6
#define NE_NVE_BIT_POS 7

#define V_PVE_MASK  (1 << V_PVE_BIT_POS)
#define V_NVE_MASK  (1 << V_NVE_BIT_POS)
#define H_PVE_MASK  (1 << H_PVE_BIT_POS)
#define H_NVE_MASK  (1 << H_NVE_BIT_POS)
#define NW_PVE_MASK (1 << NW_PVE_BIT_POS)
#define NW_NVE_MASK (1 << NW_NVE_BIT_POS)
#define NE_PVE_MASK (1 << NE_PVE_BIT_POS)
#define NE_NVE_MASK (1 << NE_NVE_BIT_POS)

/* sign coding context */
#define SC_BITS 8
#define SC_MASK (~((-1) << SC_BITS))

#endif
