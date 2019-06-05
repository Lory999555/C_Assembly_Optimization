%include "sseutils.nasm"

section .data
	

x		equ		8

y   equ     12

distance		equ		16

d		equ		20

dim     equ     4


section .bss
    
section .text

global dist32A

	dist32A:
        start

        mov eax,[ebp+x]
        mov ebx,[ebp+y]

        movaps xmm1,[eax] ;x[0]
        subps xmm1,[ebx]  ;x[0]-y[0]
        mulps xmm1,xmm1     ;(..)^2
        ;printregps xmm1
        mov edi,[ebp+d] ;d
        ;mov ecx,dim     ;4
        sub edi,16     ;d-16
        xorps xmm2, xmm2
        mov esi,4       ;i=4
    ciclo:
        cmp esi,edi     ;(j>=d-16)?
        jg resto
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

        movaps xmm0,[eax+4*esi]
        subps xmm0,[ebx+4*esi]
        mulps xmm0,xmm0
        ;printregps xmm0
        addps xmm1,xmm0
        ;printregps xmm1
        add esi,4

        movaps xmm0,[eax+4*esi]
        subps xmm0,[ebx+4*esi]
        mulps xmm0,xmm0
        ;printregps xmm0
        addps xmm1,xmm0
        ;printregps xmm1
        add esi,4

        jmp ciclo
    resto:
        mov edi, [ebp+d]
        sub edi, 4
    ciclo2:
        cmp esi, edi
        jg  fine
        movaps xmm0,[eax+4*esi] ;sommo gli ultimi elementi rimanenti
        subps xmm0,[ebx+4*esi]
        mulps xmm0,xmm0
        addps xmm2, xmm0
        add esi, 4
        jmp ciclo2

    fine:
        haddps xmm1,xmm2        ;merge di tutte le somme
        haddps xmm1,xmm1        ;|
        haddps xmm1,xmm1        ;|
        ;printregps xmm1        

        mov eax,[ebp+distance]
        movss [eax],xmm1        ;carico il nuovo valore di distance
        stop


global dist32U

	dist32U:
         start

        mov eax,[ebp+x]
        mov ebx,[ebp+y]

        movups xmm1,[eax] ;x[0]
        movups xmm7, [ebx]
        subps xmm1,xmm7  ;x[0]-y[0]
        mulps xmm1,xmm1     ;(..)^2
        ;printregps xmm1
        mov edi,[ebp+d] ;d
        ;mov ecx,dim     ;4
        sub edi,16     ;d-16
        xorps xmm2, xmm2
        mov esi,4       ;i=4
    cicloU:
        cmp esi,edi     ;(j>=d-16)?
        jg restoU
        movups xmm0,[eax+4*esi] ;x[i]
        movups xmm7, [ebx+4*esi]
        subps xmm0,xmm7  ;x[i]-y[i]
        mulps xmm0,xmm0         ;(..)^2
        ;printregps xmm0
        addps xmm1,xmm0         ;distance+=(..)^2
        ;printregps xmm1
        add esi,4               ;avanzo di indice

        movups xmm0,[eax+4*esi] ;x[i]
        movups xmm7, [ebx+4*esi]
        subps xmm0,xmm7  ;x[i]-y[i]
        mulps xmm0,xmm0         ;(..)^2
        ;printregps xmm0
        addps xmm1,xmm0         ;distance+=(..)^2
        ;printregps xmm1
        add esi,4               ;avanzo di indice

        movups xmm0,[eax+4*esi] ;x[i]
        movups xmm7, [ebx+4*esi]
        subps xmm0,xmm7  ;x[i]-y[i]
        mulps xmm0,xmm0         ;(..)^2
        ;printregps xmm0
        addps xmm1,xmm0         ;distance+=(..)^2
        ;printregps xmm1
        add esi,4               ;avanzo di indice

        movups xmm0,[eax+4*esi] ;x[i]
        movups xmm7, [ebx+4*esi]
        subps xmm0,xmm7  ;x[i]-y[i]
        mulps xmm0,xmm0         ;(..)^2
        ;printregps xmm0
        addps xmm1,xmm0         ;distance+=(..)^2
        ;printregps xmm1
        add esi,4               ;avanzo di indice

        jmp cicloU
    restoU:
        mov edi, [ebp+d]
        sub edi, 4
    ciclo2U:
        cmp esi, edi
        jg  resto2U
        movups xmm0,[eax+4*esi] ;sommo gli ultimi elementi rimanenti
        movups xmm7, [ebx+4*esi]
        subps xmm0,xmm7
        mulps xmm0,xmm0
        addps xmm2, xmm0
        add esi, 4
        jmp ciclo2U

    resto2U:
        mov edi, [ebp+d]
    ciclo3U:
        cmp esi, edi
        je  fineU
        movss xmm0,[eax+4*esi] ;sommo gli ultimi elementi rimanenti
        subss xmm0,[ebx+4*esi]
        mulss xmm0,xmm0
        addss xmm2, xmm0
        add esi, 1
        jmp ciclo3U
    fineU:
        haddps xmm1,xmm2        ;merge di tutte le somme
        haddps xmm1,xmm1        ;|
        haddps xmm1,xmm1        ;|
        ;printregps xmm1        

        mov eax,[ebp+distance]
        movss [eax],xmm1        ;carico il nuovo valore di distance
        stop

