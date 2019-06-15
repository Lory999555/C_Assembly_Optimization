%include "sseutils.nasm"

section .data
	
	align 16
    data		equ		8

    align 16
    c   equ     12

    align 16
    distance		equ		16

    i		equ		20
    n equ 24
    d equ 28
    h equ 32

    dim equ 4

    align 16
section .bss
    
section .text

global rowDistance32

	rowDistance32:

        start
        mov ecx,[ebp+c]
        mov eax,[ebp+data]
        mov edi,[ebp+d]

        mov edx,[ebp+h]                     ;h
        imul edx,dim
        add edx,eax                         ;data[h+0*n]

        mov ebx,[ebp+i]                     ;prendo i
        imul ebx,dim
        imul ebx,edi                        ;i*d
        add ecx,ebx                         ;c[i*d+0]

        movaps xmm1,[edx]                   ;data[h+0*n]
        printregps xmm1
        subps xmm1,[ecx]                    ;data[h+0*n]-c[i*d+0]
        printregps xmm1
        mulps xmm1,xmm1                     ;(..)^2
        printregps xmm1
        
        
        sub edi,dim                          ;d-4

        mov esi,4                            ;j=4
    ciclo:
        cmp esi,edi
        jge fine

        mov eax,[ebp+n]                     ;n
        imul eax,esi,4                      ;j*n
        movaps xmm0,[edx+eax]               ;data[h+j*n]
        ;printregps xmm0
        subps xmm0,[ecx+4*esi]              ;data[h+j*n] - c[i*d+j]
        ;printregps xmm0
        mulps xmm0,xmm0                     ;(..)^2
        ;printregps xmm0
        addps xmm1,xmm0                     ;distance+=(..)^2
        ;printregps xmm1
        add esi,4                           ;avanzo di indice
        
        mov eax,[ebp+n]
        imul eax,esi,4
        movaps xmm0,[edx+4*esi]
        ;printregps xmm0
        subps xmm0,[ecx+4*esi]
        ;printregps xmm0
        mulps xmm0,xmm0
        ;printregps xmm0
        addps xmm1,xmm0                     ;distance+=(..)^2
        ;printregps xmm1
        add esi,4

        jmp ciclo
    fine:
        mov eax,[ebp+n]                     ;n
        imul eax,esi,4                      ;j*n
        movaps xmm0,[edx+4*esi] ;sommo gli ultimi elementi rimanenti
        subps xmm0,[ecx+4*esi]
        mulps xmm0,xmm0
        haddps xmm1,xmm0        ;merge di tutte le somme
        haddps xmm1,xmm1        ;|
        haddps xmm1,xmm1        ;|

        mov eax,[ebp+distance]
        movss  [eax],xmm1

        stop