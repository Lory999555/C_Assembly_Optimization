%include "sseutils.nasm"

section .data


x		equ		8


y   equ     12


distance		equ		16

d		equ		20

dim equ 4


section .bss

section .text

global dist32

dist32:
    start

    mov eax,[ebp+x]
    mov ebx,[ebp+y]

    movaps xmm1,[eax] ;x[0]
    subps xmm1,[ebx]  ;x[0]-y[0]
    mulps xmm1,xmm1     ;(..)^2
    ;printregps xmm1
    mov edi,[ebp+d] ;d
    mov ecx,dim     ;4
    sub edi,ecx     ;d-4

    mov esi,4       ;i=4
ciclo:
    cmp esi,edi     ;(j>=d-4)?
    jge fine
    movaps xmm0,[eax+4*esi] ;x[i]
    subps xmm0,[ebx+4*esi]  ;x[i]-y[i]
    mulps xmm0,xmm0         ;(..)^2
    ;printregps xmm0
    addps xmm1,xmm0         ;distance+=(..)^2
    ;printregps xmm1
    add esi,4               ;avanzo di indice

    movaps xmm0,[eax+4*esi]
    subps xmm0,[ebx+4*esi]
    mulps xmm0,xmm0
    ;printregps xmm0
    addps xmm1,xmm0
    ;printregps xmm1
    add esi,4
    jmp ciclo
fine:
    movaps xmm0,[eax+4*esi] ;sommo gli ultimi elementi rimanenti
    subps xmm0,[ebx+4*esi]
    mulps xmm0,xmm0
    haddps xmm1,xmm0        ;merge di tutte le somme
    haddps xmm1,xmm1        ;|
    haddps xmm1,xmm1        ;|
    ;printregps xmm1        

    mov eax,[ebp+distance]
    movss [eax],xmm1        ;carico il nuovo valore di distance
    stop