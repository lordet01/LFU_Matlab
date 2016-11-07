#ifndef _GOLAY_H_
#define _GOLAY_H_


#include <stdio.h>


#define X22             0x00400000   /* vector representation of X^{22} */
#define X11             0x00000800   /* vector representation of X^{11} */
#define MASK12          0xfffff800   /* auxiliary vector for testing */
#define GENPOL          0x00000c75   /* generator polinomial, g(x) */

#ifdef __cplusplus
extern "C" {
#endif

__declspec(dllexport) long encode_golay(long data);
__declspec(dllexport)long decode_golay(long codeword);

#ifdef __cplusplus
}
#endif

#endif