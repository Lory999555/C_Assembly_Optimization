%include "sseutils.nasm"

section .data
	
	data equ 8
	c equ 12
	c1 equ 16
	labels equ 20
   	counts equ 24
    error equ 28
    n equ 32
    k equ 36
    d equ 40
    maxFloat dq 0x7f7fffff                ;largest float number

section .bss
	
section .text

global closestCluster32

	closestCluster32:

        start
        
        movss xmm0,[ebp+data]
        movss xmm1,[ebp+c]
        mov eax,[ebp+k]
        mov edi,[ebp+d]
        ;mov edx,[ebp+labels]
        
        mov esi,0
    forh:
        mov ebx,0
        movaps xmm4,maxFloat

    fori:
        xor xmm2,xmm2
        mov ecx,0
    forj:
        mov edx,edi
        imul edx,esi,4
        addps edx,xmm0
        movaps xmm3,[edx+ecx]

        mov edx,ebx
        imul edx,edi,4
        addps edx,xmm1
        subps xmm3,[edx+ecx]

        mulps xmm3,xmm3
        addps xmm2,xmm3
        add ecx,4
        cmp ecx,edi
        jl forj

        cmp xmm2,xmm4
        jge endif
        mov [edx+esi],ebx
        mov xmm4,xmm2
    endif:
        inc ebx
        cmp ebx,eax
        jl fori

        movss xmm1,[ebp+c1]
        mov ecx,0
        movss xmm2,[ebp+labels+esi]
    forj2:
        mov edx,edi
        mulps edx,xmm2
        addps edx,xmm1
        movaps xmm3,[edx+ecx]

        mov edx,edi
        mul edx,esi
        addps edx,xmm0
        movaps xmm4,[edx+ecx]
        addps xmm3,xmm4
        add ecx,4
        cmp ecx,edi
        jl forj2

        movss xmm3,[ebp+counts]
        inc [ebp+counts+labels+esi]
        mov ebx,[ebp+error]
        add ebx,xmm4

        inc esi
        mov ebx,[ebp+n]
        cmp esi,ebx
        jl forh


        stop