%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

align 16
list		equ		8

value       equ		12

n			equ		16


dim     equ     4
p		equ		4

section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global	assignValue

assignValue:

		push		ebp							; salva il Base Pointer
		mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
		push		ebx							; salva i registri da preservare
		push		esi
		push		edi

		mov			ebx,[ebp+n]			;n
		imul		ebx,dim				;n*4
   
        movss       xmm0,[ebp+value]          ;max_float
      
        shufps      xmm0,xmm0,0
       	;printregps  xmm0
		mov 		eax,0				;i

fori:

        mov         edi,[ebp+list]      ;list

        movaps      [eax+edi],xmm0      ;list + i*4 <- value

		add			eax,dim*p			;i + p
		cmp			eax,ebx				; i < n
		jb			fori

		pop	        edi									; ripristina i registri da preservare
		pop	        esi
		pop	        ebx
		mov	        esp, ebp							; ripristina lo Stack Pointer
		pop	        ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante