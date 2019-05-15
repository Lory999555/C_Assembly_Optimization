%include "sseutils.nasm"

section .data
	
	data equ 8
	c equ 12
	h equ 16
   	ddd equ 20
    i equ 24

    align 16
section .bss
	
section .text

global distance32

	distance32:

        push		ebp							; salva il Base Pointer
		mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
		push		ebx							; salva i registri da preservare
		push		esi
		push		edi

        mov eax,[ebp+data]
        mov ecx,[ebp+c]
        mov edi,[ebp+ddd]

        mov esi,0
        xorps xmm5,xmm5
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
        

        pop edi
        pop esi
        pop ebx
		mov	esp, ebp
		pop	ebp
		ret