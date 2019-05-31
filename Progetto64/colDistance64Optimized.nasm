; --------------------------------------------------
; Esempi di utilizzo delle istruzioni SSE
; --------------------------------------------------
; F. Angiulli
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

n			equ		16

dim		equ		4
p		equ		8
UNROLL		equ		4
BLOCKSIZE	equ		32



section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global	colDistance64OptimizedA
global	colDistance64OptimizedU

colDistance64OptimizedA:

	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali


	mov			r15,[rbp+n]		;n
	
	;mov		rcx, [rbp+starti]	; i = starti
	imul	rcx, dim		; i= i*4

	;mov		r8, [rbp+startj]		; j = startj
	imul	r8,dim		;j=j*4

	mov 	r14,r9		;d
	imul	r14,dim		;d*4

	;movaps		xmm0, [C+rcx+r10]	; c0 = C[i..i+p-1][j]  = C[dim*(i+j*n)..dim*(i+j*n+p-1)]		
	mov		r13, 0			; k = 0

	vxorps ymm4,ymm4		; xmm2=0;
	vxorps ymm5,ymm5		; xmm2=0;
	vxorps ymm6,ymm6		; xmm2=0;
	vxorps ymm7,ymm7		; xmm2=0;
	;vxorps xmm7,xmm7

fork1:		
	
	mov 		r12,rsi		;centroids
	mov 		r10,r9		;d
	imul		r10,r8		; 4*j*d
	add			r12, r10		; centroids + 4*j*d
	vmovss		xmm4, [r13+r12]	; C[j][k] = C[4*k+4*j*d]
	
	mov 		r12,rdi		;dataset
	mov			r10,r15		;n
	imul		r10, r13		; 4*k*n
	add			r12, r10		;dataset + 4*k*n

	vmovaps		ymm0, [rcx+r12]	; DS[i..i+p-1][k] = DS[4*i+4*k*n..4*(i+*k*n+p-1)]
	vmovaps		ymm1, [rcx+r12+32]	;
	vmovaps		ymm2, [rcx+r12+64]	;
	vmovaps		ymm3, [rcx+r12+96]	;




	


	vshufps		ymm4, ymm4, 0
	;printregps  xmm1
	vsubps		ymm0, ymm3
	vsubps		ymm1, ymm3
	vsubps		ymm2, ymm3
	vsubps		ymm3, ymm3

	vmulps		ymm0, ymm0		;tmp[i..i+p-1] rispetto al j-r12mo centroide
	vmulps		ymm1, ymm1
	vmulps		ymm2, ymm2
	vmulps		ymm3, ymm3
	;sqrtps		xmm0, xmm0

	vaddps		ymm4, ymm0		; distance[i..i+p-1] += tmp
	vaddps		ymm5, ymm1
	vaddps		ymm6, ymm2
	vaddps		ymm7, ymm3

	;movaps		xmm0, [rcx+r12+48]	;
	;vsubps		xmm0, xmm3
	;vmulps		xmm0, xmm0
	;vaddps		xmm7, xmm0

	add			r13, dim		; k++
	;vaddps   	xmm7,[inizio]
	cmp			r13, r14		; (k < dimension) ?
	jb			fork1
	
	;mov			r12,[rbp+distance]
	vmovaps		[rdx], ymm4	; dist[i..i+p-1]  = distance[i..i+p-1] 
	vmovaps		[rdx+32],ymm5
	vmovaps		[rdx+64],ymm6
	vmovaps		[rdx+96],ymm7


	;printregps		xmm2
	popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp				; ripristina il Base Pointer
	ret						; torna alla funzione C chiamante



colDistance64OptimizedU:

	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

		mov			r15,[rbp+n]		;n

	;mov		rcx, [rbp+starti]	; i = starti
	imul	rcx, dim		; i= i*4

	;mov		r8, [rbp+startj]		; j = startj
	imul	r8,dim		;j=j*4

	mov 	r14,r9		;d
	imul	r14,dim		;d*4

	;movaps		xmm0, [C+rcx+r10]	; c0 = C[i..i+p-1][j]  = C[dim*(i+j*n)..dim*(i+j*n+p-1)]		
	mov		r13, 0			; k = 0

	vxorps ymm4,ymm4		; xmm2=0;
	vxorps ymm5,ymm5		; xmm2=0;
	vxorps ymm6,ymm6		; xmm2=0;
	vxorps ymm7,ymm7		; xmm2=0;
	;vxorps xmm7,xmm7

forkk1:		
	;printregps xmm7
	mov 		r12,rsi		;centroids
	mov 		r10,r9		;d
	imul		r10,r8		; 4*j*d
	add			r12, r10		; centroids + 4*j*d
	vmovss		xmm4, [r13+r12]	; C[j][k] = C[4*k+4*j*d]
	
	mov 		r12,rdi		;dataset
	mov			r10,r15		;n
	imul		r10, r13		; 4*k*n
	add			r12, r10		;dataset + 4*k*n

	vmovups		ymm0, [rcx+r12]	; DS[i..i+p-1][k] = DS[4*i+4*k*n..4*(i+*k*n+p-1)]
	vmovups		ymm1, [rcx+r12+32]	;
	vmovups		ymm2, [rcx+r12+64]	;
	vmovups		ymm3, [rcx+r12+96]	;




	


	vshufps		ymm4, ymm4, 0
	;printregps  xmm1
	vsubps		ymm0, ymm4
	vsubps		ymm1, ymm4
	vsubps		ymm2, ymm4
	vsubps		ymm3, ymm4

	vmulps		ymm0, ymm0		;tmp[i..i+p-1] rispetto al j-r12mo centroide
	vmulps		ymm1, ymm1
	vmulps		ymm2, ymm2
	vmulps		ymm3, ymm3
	;sqrtps		xmm0, xmm0

	vaddps		ymm4, ymm0		; distance[i..i+p-1] += tmp
	vaddps		ymm5, ymm1
	vaddps		ymm6, ymm2
	vaddps		ymm7, ymm3

	;movaps		xmm0, [rcx+r12+48]	;
	;vsubps		xmm0, xmm3
	;vmulps		xmm0, xmm0
	;vaddps		xmm7, xmm0

	add			r13, dim		; k++
	;vaddps   	xmm7,[inizio]
	cmp			r13, r14		; (k < dimension) ?
	jb			forkk1
	
	;mov			r12,[rbp+distance]
	vmovaps		[rdx], ymm4	; dist[i..i+p-1]  = distance[i..i+p-1] 
	vmovaps		[rdx+32],ymm5
	vmovaps		[rdx+64],ymm6
	vmovaps		[rdx+96],ymm7

	;printregps		xmm2
	popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp				; ripristina il Base Pointer
	ret						; torna alla funzione C chiamante
