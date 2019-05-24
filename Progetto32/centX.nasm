%include "sseutils.nasm"

section .data
	
	align 16
    cent	equ		8
    align 16
    xx      equ     12
    k       equ     16
    dd       equ     20
    align 16
    tmp     equ     24  ;risultato di dist
    park    equ     28  ;risultato di centX
    dim1     equ     4

    section .bss
    section .text

    global centX

        centX:
            start
            mov     eax, [ebp+xx]   ;x            
            mov     ebx, [ebp+cent] ;cent
            mov     ecx, [ebp+tmp]  ;tmp
            movss   xmm5, 0.0       ; dis=0
            push    eax
            push    ebx
            push    xmm5
            call    dist32
            pop     xmm5
            pop     ebx
            pop     eax

            

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

align 16
x		    equ		20
align 16
y           equ     16
align 16
distance	equ		12
d		    equ		8
dim         equ     4

dist32:
        push    ebp
        mov     ebp,esp

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

        movss xmm5,[ebp+distance]
        movss [xmm5],xmm1        ;carico il nuovo valore di distance
        pop     ebp
        ret

