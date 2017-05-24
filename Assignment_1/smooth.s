.data
N_COEFFS: .word  3;
coeff:	.double 0.5, 1.0, 0.5;
N_SAMPLES: .word 5;
sample:	.double 1.0, 2.0, 1.0, 2.0, 1.0;
result:	.double 0.0, 0.0, 0.0, 0.0, 0.0;

.text
start: 	lwu  	R1, N_COEFFS(R0);               R1 = N_COEFFS
		lwu  	R2, N_SAMPLES(R0);				R2 = N_SAMPLES
		daddu 	R5, R0, R0;                     R5 = j = 0
		slt  	R3, R2, R1;						R3 = N_SAMPLES < N_COEFFS
		daddu   R4, R0, R0;						R4 = i = 0
		bnez 	R3, out;					    if (R6) goto endmain
		mtc1  	R0, F1;							F1 = 0
		jal  	func;							call subroutine
out:    nop
        halt

;
; parameter passed in sample,coeff,result,N_SAMPLES, return value in result
; F11 = result[i]
;
func:	mtc1  	R0, F0;							F0 = 0
sumnorm: dsll   R7, R5, 3;						R7 = j * 8
		l.d     F3, coeff(R7);					F3 = coeff[j]
	    c.lt.d 	F3, F1;							if(F1=0>F3) FP flag==TRUE
	    bc1f  	add;							if(F1=0<F3) goto add
	    daddui  R5, R5, 1; 						j++
	    sub.d   F3, F1, F3;						F3=0-F3
add:	add.d   F0, F0, F3;						norm += abs(coeff[j])
		slt 	R7, R5, R1; 					R7 = j<N_COEFFS
		bnez 	R7,	sumnorm; 					if (R7) goto sumnorm
		
		daddui  R4, R4, 1;						i = 1
		daddu   R5, R0, R0;                     R5 = j = 0
		l.d     F2, sample(R0);					F2 = sample[0]
		s.d     F2, result(R0);					result[0] = F2
		
loop:	daddu   R6, R4, R5;                     R6 = i+j
		daddui  R7, R6, -1;					    R7 = i+j-1
		dsll    R8, R5, 3;						R8 = j*8
		dsll    R7, R7, 3;						R7 = (i+j-1)*8
		l.d     F2, sample(R7);					F2 = sample[i+j-1]
		l.d     F3, coeff(R8);					F3 = coeff[j]
		daddui  R5, R5, 1;						R5 = j+1
		dsll    R7, R6, 3;						R7 = (i+j)*8
		dsll    R8, R5, 3;						R8 = (j+1)*8
		l.d     F5, sample(R7);					F5 = sample[i+j]
		l.d     F6, coeff(R8);					F6 = coeff[j+1]
		daddui  R5, R5, 1;						R5 = j+2
		daddui  R7, R6, 1;						R7 = i+j+1
		dsll    R8, R5, 3;						R8 = (j+2) * 8
		dsll    R7, R7, 3;						R7 = (i+j+1) * 8
		l.d     F8, sample(R7);					F8 = sample[i+j+1]
		l.d     F9, coeff(R8);					F9 = coeff[j+2]
		mul.d   F4, F2, F3;						F4 = sample[i+j-1]*coeff[j]
		mul.d   F7, F5, F6;						F7 = sample[i+j] * coeff[j+1]
		mul.d   F10, F8, F9;					F10 = sample[i+j+1] * coeff[j+2]
		add.d   F11, F7, F4;					F11 = sample[i+j] * coeff[j+1] + sample[i+j-1] * coeff[j]
		add.d   F11, F11, F10;					F11 = F11 + sample[i+j+1] * coeff[j+2]
		slt 	R7, R5, R1; 					R7 = (j+k-1)<N_COEFFS
		daddui  R5, R5, 1; 						j+=3
		bnez 	R7,	loop; 						if (R7) goto loop

		slt 	R7, R5, R1; 					R7 = j<N_COEFFS
		beqz 	R7,	endrestj; 					if (!R7) goto endrestj

restj:	dsll    R8, R5, 3;						R8 = j*8
		daddu   R6, R4, R5;                     R6 = i+j
		daddui  R7, R6, -1;					    R7 = i+j-1
		dsll    R7, R7, 3;						R7 = (i+j-1) * 8
		l.d     F3, coeff(R8);					F3 = coeff[j]
		l.d     F2, sample(R7);					F2 = sample[i+j-1]
		mul.d   F4, F2, F3;						F4 = sample[i+j-1] * coeff[j]
		daddui  R5, R5, 1; 						j++
		slt 	R7, R5, R1; 					R7 = j<N_COEFFS
		add.d   F11, F11, F4;					F11 = F11 + sample[i+j-1] * coeff[j]
		bnez 	R7,	restj; 						if (R7) goto restj

endrestj: div.d  F11, F11, F0;					F11 = F9/F0
		dsll    R7, R4, 3;						R7 = i*8
		daddui  R4, R4, 1; 						i++
		daddui  R6, R2, -1;						R6 = N_SAMPLES-1
		slt 	R11, R4, R6; 					R11 = i<(N_SAMPLES-1)
		daddu   R5, R0, R0;                     j = 0
		bnez 	R11, loop; 						if (R11) goto loop
		s.d     F11, result(R7);				result[i] = F11

endfor1: dsll    R6, R6, 3;						R6 = (N_SAMPLES-1)*8
		l.d     F2, sample(R6);					F2 = sample[i-1]
		s.d     F2, result(R6);					result[i-1] = F2
		jr 		R31