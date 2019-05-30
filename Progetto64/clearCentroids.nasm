%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati


counts		equ		8


centroids_1		equ		12

k		equ		16
dimension		equ		20

dim		equ		4
p		equ		8
UNROLL		equ		4
BLOCKSIZE	equ		32


zero:		dd		0.0

align 16
zero_v:     dd      0.0,0.0,0.0,0.0


;align 16
;uno:		dd		1.0, 1.0, 1.0, 1.0

;align 16
;due:		dd		0.00001, 0.00001, 0.00001, 0.00001

section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global	clearCentroids

clearCentroids:

	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	mov     eax,[ebp+k]     ;k
	imul	eax,dim			;k*4

	mov 	edx,[ebp+dimension]		;d
		
	;xorps 	xmm7,xmm7
	;xorps 	xmm6,xmm6
	movss   xmm0,[zero]
	movaps   xmm1,[zero_v]


	mov		ebx, 0			; i = 0
fori:		
	mov     esi,[ebp+counts]
	movss   [ebx+esi],xmm0
	mov		ecx, 0			; j = 0
forj:		

	mov 		esi,[ebp+centroids_1]		;centroids_1
	mov			edi,edx		;d
	imul		edi, ebx		; 4*i*d
	add			esi, edi		;centroids_1 + 4*i*d
	movaps      [ecx+esi],xmm1      ;centroids_1 + 4*i*d + 4*j <- 0;

	mov 	esi,edx			;d
	imul	esi,dim			;d*4

	add		ecx, dim*p		; j+=p
	cmp		ecx, esi		; (j < d) ?
	jb		forj
	
	add		ebx, dim		; i ++
	cmp		ebx, eax		; (i < k) ?
	jb		fori

	popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp				; ripristina il Base Pointer
	ret						; torna alla funzione C chiamante
