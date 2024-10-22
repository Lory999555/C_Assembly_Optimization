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

global	colDistance64Sing

colDistance64Sing:

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

	vxorps 	ymm4,ymm4


fork4:		

	mov 		r12,rdi		;dataset
	mov			r10,r15		;n
	imul		r10, r13		; 4*k*n
	add			r12, r10		;dataset + 4*k*n
	vmovss		xmm0, [rcx+r12]	; DS[i..i+p-1][k] = DS[4*i+4*k*n..4*(i+*k*n+p-1)]


	;printregps  xmm2
	mov 		r12,rsi		;centroids
	mov 		r10,r9		;d
	imul		r10,r8		; 4*j*d
	add			r12, r10		; centroids + 4*j*d
	vmovss		xmm3, [r13+r12]	; C[j][k] = C[4*k+4*j*d]
	


	vsubss		xmm0, xmm3
	vmulss		xmm0, xmm0		;tmp[i..i+p-1] rispetto al j-r12mo centroide
	vaddss		xmm4, xmm0		; distance[i..i+p-1] += tmp

	add			r13, dim		; k++
	;addps   	xmm7,[inizio]
	cmp			r13, r14		; (k < dimension) ?
	jb			fork4
	
	vmovss		[rdx], xmm4


	popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp				; ripristina il Base Pointer
	ret						; torna alla funzione C chiamante
