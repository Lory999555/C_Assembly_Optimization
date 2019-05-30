%include "sseutils64.nasm"

section .data


c		equ		8


x   equ     12


distance		equ		16

i		equ		20
j		equ		24
k equ 28
subb equ 32

dim equ 4


section .bss

section .text

global rowDistance32Adc

rowDistance32Adc:

    start
    mov ecx,[ebp+c]
    mov ebx,[ebp+k]
    mov edi,[ebp+subb]

    mov edx,[ebp+j]                     ;prendo j
    imul edx,dim   
    imul edx,edi                        ;j*sub
    imul edx,ebx                        ;j*k*subb
    ;printreg edx

    mov ebx,[ebp+i]                     ;prendo i
    imul ebx,dim
    imul ebx,edi                        ;i*subb
    add ebx,edx                         ;j*k*subb+i*subb
    add ebx,ecx                         ;c[j*k*subb+i*subb]

    mov eax,[ebp+x]                     ;prendo x
    
    
    ;xorps xmm5,xmm5                     ;distance=0 

    movaps xmm1,[eax]                   ;uj_x[0]
    subps xmm1,[ebx]                    ;uj_x[0]-c[0*k*subb+i*subb]
    mulps xmm1,xmm1                     ;(..)^2
    ;printregps xmm1
    mov edi,[ebp+subb]                     ;subb
    mov ecx,dim
    sub edi,ecx                          ;subb-4

    mov esi,4                            ;z=4
ciclo:
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;movaps xmm0,[eax+4*esi]             ;uj_x[z] (prendo 4 elementi la volta)
    ;printregps xmm0

    ;movaps xmm1,[edx+4*esi]             ;c[j*k*subb+i*subb+z]  (prendo 4 elementi la volta)
    
    ;printregps xmm1

    ;subps xmm0,xmm1                     ;uj_x[z] - c[j*k*subb+i*subb+z]

    ;mulps xmm0,xmm0                     ;(...)^2
    
    ;xorps xmm2,xmm2                     ;tmp_dist=0
    ;movss xmm2,xmm5                     ;carico distance
    ;printregps xmm2
    ;haddps xmm2,xmm0                    ;sommo la nuova distanza con la vecchia distanza

    ;haddps xmm2,xmm2

    ;haddps xmm2,xmm2                    ;avrò un vettore con 4 valori uguali

    ;movss xmm5,xmm2                     ;carico uno dei valori di xmm2 in xmm5

    ;add esi,4
    ;cmp esi,edi
    ;jl ciclo
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    cmp esi,edi
    jge fine

    movaps xmm0,[eax+4*esi]             ;uj_x[z]
    ;printregps xmm0
    subps xmm0,[ebx+4*esi]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps xmm0
    mulps xmm0,xmm0                     ;(..)^2
    ;printregps xmm0
    addps xmm1,xmm0                     ;distance+=(..)^2
    ;printregps xmm1
    add esi,4                           ;avanzo di indice

    movaps xmm0,[eax+4*esi]
    ;printregps xmm0
    subps xmm0,[ebx+4*esi]
    ;printregps xmm0
    mulps xmm0,xmm0
    ;printregps xmm0
    addps xmm1,xmm0                     ;distance+=(..)^2
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

    mov eax,[ebp+distance]
    movss  [eax],xmm1

    stop