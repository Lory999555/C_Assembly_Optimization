%include "sseutils.nasm"

section .data
	
	align 16
    x		equ		8

    align 16
    y   equ     12

    align 16
    distance		equ		16

    d		equ		20
    i equ 24
    dim equ 4

    align 16
section .bss
    
section .text

global wncDist32

	wncDist32:
        start

        mov eax,[ebp+x]
        mov ebx,[ebp+y]
        mov edx,[ebp+i]     ;i
        mov edi,[ebp+d] ;d
        imul edx,edi ;i*d
        add ebx,edx         ;y[i*d]

        movaps xmm1,[eax] ;x[0]
        ;printregps xmm1
        subps xmm1,[ebx]  ;x[0]-y[i*d]
        ;printregps xmm1
        mulps xmm1,xmm1     ;(..)^2
        ;printregps xmm1

        mov ecx,dim     ;4
        sub edi,ecx     ;d-4

        mov esi,4       ;j=4
    ciclo:
        cmp esi,edi     ;(j>=d-4)?
        jge fine
        movaps xmm0,[eax+4*esi] ;x[j]
        ;printregps xmm0
        subps xmm0,[ebx+4*esi]  ;x[j]-y[i*d+j]
        ;printregps xmm0
        mulps xmm0,xmm0         ;(..)^2
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