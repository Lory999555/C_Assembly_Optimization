%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati



dim		equ		4
p		equ		8
UNROLL		equ		4
BLOCKSIZE	equ		32




;align 16
;uno:		dd		1.0, 1.0, 1.0, 1.0

;align 16
;due:		dd		0.00001, 0.00001, 0.00001, 0.00001

section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global	updateCentroid

updateCentroid:

	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	;mov     rcx,[rbp+k]     ;k
	imul	rcx,dim			;k*4

	;mov 	r8,[rbp+dimension]		;d
		
	vxorps 	ymm15,ymm15
	;xorps 	xmm6,xmm6
	


	mov		rbx, 0			; i = 0
fori:		
	cmp		r8,p			;d < 8
	jb		forj1
	mov		r13, 0			; j = 0
forj:		
	;attenzione ad r11 che viene usato dopo quindi se lo si modifica potrebbe non funzionare
	;mov         r11,[rbp+counts]    ;counts

	vmovss		xmm1,[rdx+rbx]		;counts + i*4
	vcomiss     xmm1,xmm15          ; counts + i*4 != 0
	jz         else

	;printregps	xmm7
	;addps		xmm7,[uno]

	;movss       xmm1,[r11+rbx]      ;espando il counts[i] su xmm1
	;movss       xmm1,[r11+rbx]
	vbroadcastss ymm1,[rdx+rbx]   ;counts[i] ripetuto 8 volte
	;printregps	xmm1
	
	;mov 		r11,[rbp+centroids_1]		;centroids_1
	mov 		r11,rsi		;centroids_1
	mov			r10,r8		;d
	imul		r10, rbx		; 4*i*d
	add			r11, r10		;centroids_1 + 4*i*d
	vmovups      ymm0,[r13+r11]      ;centroids_1 + 4*i*d + 4*j. prendo 8 elementi (e dimensioni)
	vdivps       ymm0,ymm1       ;c1[i*d+j..i+*d+j+p-1]/counts[i]
	;addps		xmm6,[uno]
	;printregps	xmm6

	;mov         r11,[rbp+centroids]         ;centroids
	mov         r11,rdi         ;centroids
	add         r11,r10         ;centroids + 4*i*d
	vmovups		[r13+r11], ymm0	        ; centroids + 4*i*d + 4*j <- xmm0
	jmp iter
else:
	;printregps	xmm6
	;addps		xmm6,[due]

	;mov 		r11,[rbp+centroids_1]		;centroids_1
	mov 		r11,rsi		;centroids_1
	mov			r10,r8		;d
	imul		r10, rbx		; 4*i*d
	add			r11, r10		;centroids_1 + 4*i*d
	vmovups      ymm0,[r13+r11]      ;centroids_1 + 4*i*d + 4*j. prendo 8 elementi (e dimensioni)


	mov         r11,rdi         ;centroids
	add         r11,r10         ;centroids + 4*i*d
	vmovups		[r13+r11], ymm0	        ; centroids + 4*i*d + 4*j = xmm0

iter:
	mov 	r11,r8			;d
	sub		r11,p		;d-p
	imul	r11,dim		;d*4

	add		r13, dim*p		; j+=p
	cmp		r13, r11		; (j <= d-p) ?
	jle		forj

forj1:
	vmovss		xmm1,[rdx+rbx]		;counts + i*4
	vcomiss     xmm1,xmm15          ; counts + i*4 != 0
	jz         else1

	;printregps	xmm7
	;addps		xmm7,[uno]

	;movss       xmm1,[r11+rbx]      ;espando il counts[i] su xmm1
	;movss       xmm1,[r11+rbx]
	;vbroadcastss ymm1,[rdx+rbx]   ;counts[i] ripetuto 8 volte
	;printregps	xmm1
	;mov 		r11,[rbp+centroids_1]		;centroids_1
	mov 		r11,rsi		;centroids_1
	mov			r10,r8		;d
	imul		r10, rbx		; 4*i*d
	add			r11, r10		;centroids_1 + 4*i*d
	vmovss      xmm0,[r13+r11]      ;centroids_1 + 4*i*d + 4*j. prendo 8 elementi (e dimensioni)

	vdivss       xmm0,xmm1       ;c1[i*d+j..i+*d+j+p-1]/counts[i]
	;addps		xmm6,[uno]
	;printregps	xmm6

	;mov         r11,[rbp+centroids]         ;centroids
	mov         r11,rdi         ;centroids
	add         r11,r10         ;centroids + 4*i*d
	vmovss		[r13+r11], xmm0	        ; centroids + 4*i*d + 4*j <- xmm0
	jmp iter1
else1:
	;printregps	xmm6
	;addps		xmm6,[due]

	;mov 		r11,[rbp+centroids_1]		;centroids_1
	mov 		r11,rsi		;centroids_1
	mov			r10,r8		;d
	imul		r10, rbx		; 4*i*d
	add			r11, r10		;centroids_1 + 4*i*d
	vmovss      xmm0,[r13+r11]      ;centroids_1 + 4*i*d + 4*j. prendo 8 elementi (e dimensioni)


	mov         r11,rdi         ;centroids
	add         r11,r10         ;centroids + 4*i*d
	vmovss		[r13+r11], xmm0	        ; centroids + 4*i*d + 4*j = xmm0

iter1:
	mov 	r11,r8			;d
	imul	r11,dim		;d*4

	add		r13, dim		; j+=1
	cmp		r13, r11		; (j < d) ?
	jb		forj1

	add		rbx, dim		; i ++
	cmp		rbx, rcx		; (i < k) ?
	jb		fori

	popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp				; ripristina il Base Pointer
	ret						; torna alla funzione C chiamante
