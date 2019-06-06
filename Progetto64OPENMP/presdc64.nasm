%include "sseutils64.nasm"



section .data			; Sezione contenente dati inizializzati

c1		equ		8
c2		equ		12
dim		equ		16
dist    equ     20

p		equ		8
UNROLL		equ		4
BLOCKSIZE	equ		32

section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global	rawdistance32

rawdistance32:

    push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

    mov eax,[ebp+c1]
    mov ebx,[ebp+c2]
    mov edi,[ebp+d]

    mov esi,0       ;i=0
forj:
    mov edx,[ebp+h]
    imul edx,edi
    add edx,eax
    movaps xmm0,[edx+4*esi]
    
    mov edx,[ebp+i]
    imul edx,edi
    add edx,ecx
    movaps xmm1,[edx+4*esi]

    subps xmm0,xmm1
    mulps xmm0,xmm0
    ; somma tutti i singoli
    
    movaps xmm2,xmm5
    
    haddps xmm2,xmm0

    ;printregps xmm2
    haddps xmm2,xmm2

    ;printregps xmm2
    haddps xmm2,xmm2

    ;printregps xmm2
    movaps xmm5,xmm2
    ;printregps xmm5

    add esi,4
    cmp esi,edi
    jl forj

    printregps xmm5
    extractps  eax,xmm5,00000000
    

    popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp				; ripristina il Base Pointer
	ret						; torna alla funzione C chiamante

