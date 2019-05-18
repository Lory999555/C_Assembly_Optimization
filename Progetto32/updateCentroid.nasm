%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

align 16
centroids		equ		8

align 16
centroids_1		equ		12

align 16
counts		equ		16

k		equ		20
dimension		equ		24

dim		equ		4
p		equ		4
UNROLL		equ		4
BLOCKSIZE	equ		32


align 16
zero:		dd		0.0

;align 16
;uno:		dd		1.0, 1.0, 1.0, 1.0

;align 16
;due:		dd		0.00001, 0.00001, 0.00001, 0.00001

section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global	updateCentroid

updateCentroid:

		push		ebp							; salva il Base Pointer
		mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
		push		ebx							; salva i registri da preservare
		push		esi
		push		edi

        mov     eax,[ebp+k]     ;k
		imul	eax,dim			;k*4

		mov 	edx,[ebp+dimension]		;d
		    
		;xorps 	xmm7,xmm7
		;xorps 	xmm6,xmm6


		mov		ebx, 0			; i = 0
fori:		
		mov		ecx, 0			; j = 0
forj:		
		;attenzione ad esi che viene usato dopo quindi se lo si modifica potrebbe non funzionare
        mov         esi,[ebp+counts]    ;counts
		movss		xmm1,[esi+ebx]		;counts + i*4
        comiss      xmm1,[zero]           ; counts + i*4 != 0
        jz         else

		;printregps	xmm7
		;addps		xmm7,[uno]

        ;movss       xmm1,[esi+ebx]      ;espando il counts[i] su xmm1
		;movss       xmm1,[esi+ebx]
        shufps      xmm1,xmm1,0   ;counts[i] ripetuto 4 volte
		;printregps	xmm1
        
      	mov 		esi,[ebp+centroids_1]		;centroids_1
		mov			edi,edx		;d
		imul		edi, ebx		; 4*i*d
		add			esi, edi		;centroids_1 + 4*i*d
        movaps      xmm0,[ecx+esi]      ;centroids_1 + 4*i*d + 4*j. prendo 4 elementi (e dimensioni)

        divps       xmm0,xmm1       ;c1[i*d+j..i+*d+j+p-1]/counts[i]
		;addps		xmm6,[uno]
		;printregps	xmm6

        mov         esi,[ebp+centroids]         ;centroids
        add         esi,edi         ;centroids + 4*i*d
		movaps		[ecx+esi], xmm0	        ; centroids + 4*i*d + 4*j <- xmm0
        jmp iter
else:
		;printregps	xmm6
		;addps		xmm6,[due]

        mov 		esi,[ebp+centroids_1]		;centroids_1
		mov			edi,edx		;d
		imul		edi, ebx		; 4*i*d
		add			esi, edi		;centroids_1 + 4*i*d
        movaps      xmm0,[ecx+esi]      ;centroids_1 + 4*i*d + 4*j. prendo 4 elementi (e dimensioni)


        mov         esi,[ebp+centroids]         ;centroids
        add         esi,edi         ;centroids + 4*i*d
        movaps		[ecx+esi], xmm0	        ; centroids + 4*i*d + 4*j = xmm0

iter:
		mov 	esi,edx			;d
		imul	esi,dim			;d*4

		add		ecx, dim*p		; j+=p
		cmp		ecx, esi		; (j < d) ?
		jb		forj
		
		add		ebx, dim		; i ++
		cmp		ebx, eax		; (i < k) ?
		jb		fori

		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp							; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante