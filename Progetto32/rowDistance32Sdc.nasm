%include "sseutils.nasm"

section .data
	
	align 16
    c		equ		8

    align 16
    distance		equ		12

    i		equ		16
    j		equ		20
    j_d equ 24
    k equ 28
    subb equ 32

    dim equ 4

    align 16
section .bss
    
section .text

global rowDistance32Sdc

	rowDistance32Sdc:

        start
        mov ecx,[ebp+c]
        mov eax,[ebp+k]
        mov edi,[ebp+subb]

       
        ;mov dword [distance],0                    ;distance=0
        
        mov edx,[ebp+j]                     ;prendo j
        imul edx,dim   
        imul edx,edi                        ;j*sub
        imul edx,eax                        ;j*k*subb

        mov ebx,[ebp+i]                     ;prendo i
        imul ebx,dim
        imul ebx,edi                        ;i*subb
        add ebx,edx                         ;j*k*subb+i*subb
        add ebx,ecx                         ;c[j*k*subb+i*subb]

        mov eax,[ebp+j_d]                   ;prendo j_d
        imul eax,dim
        imul eax,edi                        ;j_d*subb
        add eax,edx                         ;j*k*subb+j_d*subb
        add eax,ecx                         ;c[j*k*subb+j_d*subb]
        
        ;printreg ebx
        ;printreg eax

        xorps xmm5,xmm5                     ;distance=0 
        ;printreg edi
        mov esi,0                           ;z=0
    forj:
        
        movaps xmm0,[ebx+4*esi]             ;c[j*k*subb+i*subb+z] (prendo 4 elementi la volta)
        ;printregps xmm0

        movaps xmm1,[eax+4*esi]             ;c[j*k*subb+j_d*subb+z]  (prendo 4 elementi la volta)
        
        ;printregps xmm1

        subps xmm0,xmm1                     ;c[j*k*subb+i*subb+z] - c[j*k*subb+j_d*subb+z]

        mulps xmm0,xmm0                     ;(...)^2
        
        xorps xmm2,xmm2                     ;tmp_dist=0
        movss xmm2,xmm5                     ;carico distance
        ;printregps xmm2
        haddps xmm2,xmm0                    ;sommo la nuova distanza con la vecchia distanza

        haddps xmm2,xmm2

        haddps xmm2,xmm2                    ;avrò un vettore con 4 valori uguali

        ;printregps xmm2

        ;movss [distance],xmm2
        movss xmm5,xmm2                     ;carico uno dei valori di xmm2 in xmm5

        add esi,4
        cmp esi,edi
        jl forj

        ;movss xmm1,[distance]
        ;printregps xmm5
        mov eax,[ebp+distance]
        movss  [eax],xmm5,

        stop