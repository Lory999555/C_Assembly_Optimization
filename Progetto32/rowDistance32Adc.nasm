%include "sseutils.nasm"

section .data
	
	align 16
    c		equ		8

    align 16
    x   equ     12

    align 16
    distance		equ		16

    i		equ		20
    j		equ		24
    k equ 28
    subb equ 32

    dim equ 4

    align 16
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
        add edx,ebx                         ;j*k*subb+i*subb
        add edx,ecx                         ;c[j*k*subb+i*subb]

        mov eax,[ebp+x]                     ;prendo x
        ;printreg edx
        ;printreg eax

        mov ebx,[ebp+k]
        sub ebx,4
        xorps xmm5,xmm5                     ;distance=0 

        ;mov esi,0                           ;z=0
        mov esi,4                            ;z=4
    forj:
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        movaps xmm0,[eax+4*esi]             ;uj_x[z] (prendo 4 elementi la volta)
        ;printregps xmm0

        movaps xmm1,[edx+4*esi]             ;c[j*k*subb+i*subb+z]  (prendo 4 elementi la volta)
        
        ;printregps xmm1

        subps xmm0,xmm1                     ;c[j*k*subb+i*subb+z] - c[j*k*subb+j_d*subb+z]

        mulps xmm0,xmm0                     ;(...)^2
        
        xorps xmm2,xmm2                     ;tmp_dist=0
        movss xmm2,xmm5                     ;carico distance
        ;printregps xmm2
        haddps xmm2,xmm0                    ;sommo la nuova distanza con la vecchia distanza

        haddps xmm2,xmm2

        haddps xmm2,xmm2                    ;avrò un vettore con 4 valori uguali

        movss xmm5,xmm2                     ;carico uno dei valori di xmm2 in xmm5
    
        add esi,4
        cmp esi,edi
        jl forj
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        ;cmp esi,ebx
        ;jge fine
        ;addps xmm0,[eax+4*esi]
        ;subps xmm0,[edx+4*esi]
        ;mulps xmm0,xmm0
        ;add esi,4
        ;addps xmm0,[eax+4*esi]
        ;subps xmm0,[edx+4*esi]
        ;mulps xmm0,xmm0
        ;add esi,4

    fine:
    
        mov eax,[ebp+distance]
        movss  [eax],xmm5,

        stop