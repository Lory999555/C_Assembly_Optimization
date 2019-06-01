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


		;mov			rbx,[rbp+n]			;n
		imul		rdx,dim				;n*4
   
        ;vmovss       xmm0,[rbp+value]          ;max_float
		;vmovss       xmm0,[rbp+value]          ;max_float
      
        vbroadcastss      ymm0,[rsi]
       	;printregps  xmm0
		mov 		rax,0				;i

fori:

        ;mov         edi,[rbp+list]      ;list

        vmovaps      [rax+rdi],ymm0      ;list + i*4 <- value

		add			rax,dim*p			;i + p
		cmp			rax,rbx				; i < n
		jb			fori

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
