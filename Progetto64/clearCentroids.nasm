%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati



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

	;mov     	rax,[rbp+k]     ;k
	imul		rdx,4			;k*4

	;mov 		r12,[rbp+dimension]		;d
		
	;xorps 	xmm7,xmm7
	;xorps 	xmm6,xmm6
	vxorps  	ymm0,ymm0
	vxorps   	ymm1,ymm1


	mov			rbx, 0			; i = 0
fori:		
	;mov     	r11,[rbp+counts]
	vmovss   	[rbx+rdi],xmm0
	mov			r13, 0			; j = 0
forjclear:		

	;mov 		r11,[rbp+centroids_1]		;centroids_1
	mov			r10,rcx		;d
	imul		r10, rbx		; 4*i*d
	add			r10,rsi		;centroids_1 + 4*i*d
	vmovaps      [r13+r10],ymm1      ;centroids_1 + 4*i*d + 4*j <- 0;

	mov 	r11,rcx			;d
	imul	r11,4			;d*4

	add		r13, 32		; j+=p
	cmp		r13, r11		; (j < d) ?
	jb		forjclear
	
	add		rbx, 4		; i ++
	cmp		rbx, rdx		; (i < k) ?
	jb		fori

	popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp				; ripristina il Base Pointer
	ret						; torna alla funzione C chiamante
