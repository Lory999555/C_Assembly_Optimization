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
        movaps xmm3,[eax+4*esi+16]
        movaps xmm4,[eax+4*esi+32]
        movaps xmm5,[eax+4*esi+48]

        subps xmm0,[ebx+4*esi]  ;x[i]-y[i]
        subps xmm3,[ebx+4*esi+16]
        subps xmm4,[ebx+4*esi+32]
        subps xmm5,[ebx+4*esi+48]

        mulps xmm0,xmm0         ;(..)^2
        mulps xmm3,xmm3
        mulps xmm4,xmm4
        mulps xmm5,xmm5

        addps xmm1,xmm0         ;distance+=(..)^2
        addps xmm1,xmm3
        addps xmm1,xmm4
        addps xmm1,xmm5
        add esi,16

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
        movups xmm3,[eax+4*esi+16]
        movups xmm4,[eax+4*esi+32]
        movups xmm5,[eax+4*esi+48]

        movups xmm6,[ebx+4*esi]             ;uj_x[z]
        movups xmm7,[ebx+4*esi+16]
        

        subps xmm0,xmm6              ;uj_x[z]- c[j*k*subb+i*subb+z]
        subps xmm3,xmm7

        movups xmm6,[ebx+4*esi+32]
        movups xmm7,[ebx+4*esi+48]

        subps xmm4,xmm6
        subps xmm5,xmm7

        mulps xmm0,xmm0         ;(..)^2
        mulps xmm3,xmm3
        mulps xmm4,xmm4
        mulps xmm5,xmm5

        addps xmm1,xmm0         ;distance+=(..)^2
        addps xmm1,xmm3
        addps xmm1,xmm4
        addps xmm1,xmm5
        add esi,16

        jmp cicloU
    restoU:
        mov edi, [ebp+d]
        sub edi, 4
    ciclo2U:
        cmp esi, edi
        jg  resto2U
        movups xmm0,[eax+4*esi] ;sommo gli ultimi elementi rimanenti
        movups xmm7,[ebx+4*esi]
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

