#include <cstdint>

static const uint64_t i1 = 1;
static const uint64_t i0 = 0;

static const uint64_t primpoly[65] = {
/* 00 */                        0,
/* 01 */                       03,
/* 02 */                       07,
/* 03 */	              013,
/* 04 */	              023,
/* 05 */	              045,
/* 06 */	             0103,
/* 07 */	             0203,
/* 08 */	             0435,
/* 09 */	            01021,
/* 10 */	            02011,
/* 11 */	            04005,
/* 12 */                   010123,
/* 13 */	           020033,
/* 14 */	           040053,
/* 15 */	          0100003,
/* 16 */	          0200055,
/* 17 */	          0400011,
/* 18 */	         01000201,
/* 19 */	         02000047,
/* 20 */	         04000011,
/* 21 */	        010000005,
/* 22 */	        020000003,
/* 23 */	        040000041,
/* 24 */	       0100000033,
/* 25 */	       0200000011,
/* 26 */	       0400000107,
/* 27 */	      01000000047,
/* 28 */	      02000000011,
/* 29 */	      04000000005,
/* 30 */             010000000123,
/* 31 */	     020000000011,
/* 32 */	     040000000305,
/* 33 */	    0100000020001,
/* 34 */	    0200000000431,
/* 35 */            0400000000005,
/* 36 */	   01000000004001,
/* 37 */	   02000000000123,
/* 38 */   	   04000000000143,
/* 39 */	  010000000000021,
/* 40 */	  020000000000071,
/* 41 */	  040000000000011,
/* 42 */	 0100000000000231,
/* 43 */	 0200000000000131,
/* 44 */         0400000000000145,
/* 45 */        01000000000000033,
/* 46 */        02000000000000701,
/* 47 */        04000000000000041,
/* 48 */       010000000000001221,
/* 49 */       020000000000001001,
/* 50 */       040000000000000035,
/* 51 */      0100000000000000113,
/* 52 */      0200000000000000011,
/* 53 */      0400000000000000107,
/* 54 */     01000000000000000511,
/* 55 */     02000000000100000001,
/* 56 */     04000000000000000225,
/* 57 */    010000000000000000201,
/* 58 */    020000000000002000001,
/* 59 */    040000000000000000225,
/* 60 */   0100000000000000000003,
/* 61 */   0200000000000000000047,
/* 62 */   0400000000000000000151,
/* 63 */  01000000000000000000003,
/* 64 */   0000000000000000000033}; // 02000000000000000000033 */

static const uint64_t pw[65] = {
/* 00 */ i1,
/* 01 */ i1 << 1,
/* 02 */ i1 << 2,
/* 03 */ i1 << 3,
/* 04 */ i1 << 4,
/* 05 */ i1 << 5,
/* 06 */ i1 << 6,
/* 07 */ i1 << 7,
/* 08 */ i1 << 8,
/* 09 */ i1 << 9,
/* 10 */ i1 << 10,
/* 11 */ i1 << 11,
/* 12 */ i1 << 12,
/* 13 */ i1 << 13,
/* 14 */ i1 << 14,
/* 15 */ i1 << 15,
/* 16 */ i1 << 16,
/* 17 */ i1 << 17,
/* 18 */ i1 << 18,
/* 19 */ i1 << 19,
/* 20 */ i1 << 20,
/* 21 */ i1 << 21,
/* 22 */ i1 << 22,
/* 23 */ i1 << 23,
/* 24 */ i1 << 24,
/* 25 */ i1 << 25,
/* 26 */ i1 << 26,
/* 27 */ i1 << 27,
/* 28 */ i1 << 28,
/* 29 */ i1 << 29,
/* 30 */ i1 << 30,
/* 31 */ i1 << 31,
/* 32 */ i1 << 32,
/* 33 */ i1 << 33,
/* 34 */ i1 << 34,
/* 35 */ i1 << 35,
/* 36 */ i1 << 36,
/* 37 */ i1 << 37,
/* 38 */ i1 << 38,
/* 39 */ i1 << 39,
/* 40 */ i1 << 40,
/* 41 */ i1 << 41,
/* 42 */ i1 << 42,
/* 43 */ i1 << 43,
/* 44 */ i1 << 44,
/* 45 */ i1 << 45,
/* 46 */ i1 << 46,
/* 47 */ i1 << 47,
/* 48 */ i1 << 48,
/* 49 */ i1 << 49,
/* 50 */ i1 << 50,
/* 51 */ i1 << 51,
/* 52 */ i1 << 52,
/* 53 */ i1 << 53,
/* 54 */ i1 << 54,
/* 55 */ i1 << 55,
/* 56 */ i1 << 56,
/* 57 */ i1 << 57,
/* 58 */ i1 << 58,
/* 59 */ i1 << 59,
/* 60 */ i1 << 60,
/* 61 */ i1 << 61,
/* 62 */ i1 << 62,
/* 63 */ i1 << 63,
/* 64 */ i1 << 63}; // not used

static const uint64_t pwm1[65] = {
/* 00 */ 0,
/* 01 */ (i1 << 1) - 1,
/* 02 */ (i1 << 2) - 1,
/* 03 */ (i1 << 3) - 1,
/* 04 */ (i1 << 4) - 1,
/* 05 */ (i1 << 5) - 1,
/* 06 */ (i1 << 6) - 1,
/* 07 */ (i1 << 7) - 1,
/* 08 */ (i1 << 8) - 1,
/* 09 */ (i1 << 9) - 1,
/* 10 */ (i1 << 10) - 1,
/* 11 */ (i1 << 11) - 1,
/* 12 */ (i1 << 12) - 1,
/* 13 */ (i1 << 13) - 1,
/* 14 */ (i1 << 14) - 1,
/* 15 */ (i1 << 15) - 1,
/* 16 */ (i1 << 16) - 1,
/* 17 */ (i1 << 17) - 1,
/* 18 */ (i1 << 18) - 1,
/* 19 */ (i1 << 19) - 1,
/* 20 */ (i1 << 20) - 1,
/* 21 */ (i1 << 21) - 1,
/* 22 */ (i1 << 22) - 1,
/* 23 */ (i1 << 23) - 1,
/* 24 */ (i1 << 24) - 1,
/* 25 */ (i1 << 25) - 1,
/* 26 */ (i1 << 26) - 1,
/* 27 */ (i1 << 27) - 1,
/* 28 */ (i1 << 28) - 1,
/* 29 */ (i1 << 29) - 1,
/* 30 */ (i1 << 30) - 1,
/* 31 */ (i1 << 31) - 1,
/* 32 */ (i1 << 32) - 1,
/* 33 */ (i1 << 33) - 1,
/* 34 */ (i1 << 34) - 1,
/* 35 */ (i1 << 35) - 1,
/* 36 */ (i1 << 36) - 1,
/* 37 */ (i1 << 37) - 1,
/* 38 */ (i1 << 38) - 1,
/* 39 */ (i1 << 39) - 1,
/* 40 */ (i1 << 40) - 1,
/* 41 */ (i1 << 41) - 1,
/* 42 */ (i1 << 42) - 1,
/* 43 */ (i1 << 43) - 1,
/* 44 */ (i1 << 44) - 1,
/* 45 */ (i1 << 45) - 1,
/* 46 */ (i1 << 46) - 1,
/* 47 */ (i1 << 47) - 1,
/* 48 */ (i1 << 48) - 1,
/* 49 */ (i1 << 49) - 1,
/* 50 */ (i1 << 50) - 1,
/* 51 */ (i1 << 51) - 1,
/* 52 */ (i1 << 52) - 1,
/* 53 */ (i1 << 53) - 1,
/* 54 */ (i1 << 54) - 1,
/* 55 */ (i1 << 55) - 1,
/* 56 */ (i1 << 56) - 1,
/* 57 */ (i1 << 57) - 1,
/* 58 */ (i1 << 58) - 1,
/* 59 */ (i1 << 59) - 1,
/* 60 */ (i1 << 60) - 1,
/* 61 */ (i1 << 61) - 1,
/* 62 */ (i1 << 62) - 1,
/* 63 */ 0x7fffffffffffffff,
/* 64 */ 0xffffffffffffffff};

