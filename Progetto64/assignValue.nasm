%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati


list		equ		8

value       equ		12

n			equ		16


dim     equ     4
p		equ		8

section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global	assignValue

assignValue:

		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali


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

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
