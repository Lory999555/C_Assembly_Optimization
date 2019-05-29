%include "sseutils.nasm"

section .data
	
	align 16
    cent	equ		8
    align 16
    xx      equ     12
    k       equ     16
    ddd       equ     20
    align 16
    tmp     equ     24  ;risultato di dist
    park    equ     28  ;risultato di centX
    dis     equ     32  ;dis temporanea
    dim1     equ     4

    section .bss
    section .text

    global cent_X

        cent_X:
            start
            mov     eax, [ebp+xx]   ;x            
            mov     ebx, [ebp+cent] ;cent
            mov     ecx, [ebp+ddd]  ;d
            mov     edx, [ebp+dis]  ;dis
            mov     edi,    0
            ;printreg eax
            ;printreg ebx
            printreg edx
            ;printreg ecx
            push    eax
            push    ebx
            push    edx
            push    ecx
            call    dist32
            pop     ecx
            pop     edx
            pop     ebx
            pop     eax
            printreg edx
            printreg [edx]
            mov     edx, [ebp+dis]
            xorps   xmm2, xmm2
            movss   xmm2, [edx]

        fori: 
            mov     esi, [ebp+cent]
            mov     edx, [ebp+ddd]
            imul     edx, edi       ;i*d
            imul     edx, dim1      ;4*i*d
            add     esi, edx        ;cent + 4*i*d
            mov     ebx, esi 
            mov     ecx, [ebp+tmp]  ;tmp
            mov     edx, [ebp+ddd]
            ;mov esi, [ecx]
            printreg [ecx]
            ;printreg eax
            ;printreg ebx
            ;printreg ecx
            ;printreg edx
            push    eax
            push    ebx
            push    ecx
            push    edx
            call    dist32
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax
           ; mov esi, [ecx]
            printreg ecx
            printreg [ecx]

            mov     esi, [ebp+cent]
            mov     edx, edi
            add     edx, 1
            imul     edx, [ebp+ddd]
            imul     edx, dim1
            add     esi, edx
            mov     ebx, esi ;cent[0]
            mov     ecx, [ebp+tmp+4]
            mov     edx, [ebp+ddd]
            ;mov esi, [ecx]
            printreg [ecx]
            ;printreg eax
            ;printreg ebx
            ;printreg ecx
            ;printreg edx
            push    eax
            push    ebx
            push    ecx
            push    edx
            call    dist32
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax
            ;mov esi, [ecx]
            printreg [ecx]

            mov     esi, [ebp+cent]
            mov     edx, edi
            add     edx, 2
            imul     edx, [ebp+ddd]
            imul     edx, dim1
            add     esi, edx
            mov     ebx, esi ;cent[0]
            mov     ecx, [ebp+tmp+8]
            mov     edx, [ebp+ddd]
            ;printreg eax
            ;printreg ebx
            ;printreg ecx
            ;printreg edx
            push    eax
            push    ebx
            push    ecx
            push    edx
            call    dist32
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax

            mov     esi, [ebp+cent]
            mov     edx, edi
            add     edx, 3
            imul     edx, [ebp+ddd]
            imul     edx, dim1
            add     esi, edx
            mov     ebx, esi ;cent[0]
            mov     ecx, [ebp+tmp+12]
            mov     edx, [ebp+ddd]
            ;printreg eax
            ;printreg ebx
            ;printreg ecx
            ;printreg edx
            push    eax
            push    ebx
            push    ecx
            push    edx
            ;printreg    edi
            call    dist32
            ;printreg    edi
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax

            mov     ecx, [ebp+tmp]
            xorps   xmm1, xmm1
            movss   xmm1, [ecx]
            comiss  xmm1, xmm2
            jge      avanti
            movss   xmm2, xmm1
            mov     esi, [ebp+park]
            mov     [esi], edi
        avanti:
            mov     ecx, [ebp+tmp+4]
            movss   xmm1, [ecx]
            comiss  xmm1, xmm2
            jge      avanti2
            movss   xmm2, xmm1
            mov     esi, [ebp+park]
            mov     [esi], edi
        avanti2:
            mov     ecx, [ebp+tmp+8]
            movss   xmm1, [ecx]
            comiss  xmm1, xmm2
            jge      avanti3
            movss   xmm2, xmm1
            mov     esi, [ebp+park]
            mov     [esi], edi
        avanti3:
            mov     ecx, [ebp+tmp+12]
            movss   xmm1, [ecx]
            comiss  xmm1, xmm2
            jge      avanti4
            movss   xmm2, xmm1
            mov     esi, [ebp+park]
            mov     [esi], edi
        avanti4:
            add     edi,4 ;;;;;;;;;;;;;;;;;;farlo generale con unroll
            mov     esi, [ebp+k]
            cmp      edi, esi
            jb      fori
            stop



            

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

align 16
x		    equ		20
align 16
y           equ     16
align 16
distance	equ		12
dddd		    equ		8
dim         equ     4

dist32:
        push    ebp
        mov     ebp,esp
        push    edi
        ;push    esi
        mov eax,[ebp+distance]
        mov ebx,[ebp+dddd]
        ;printreg ebx
        ;printreg eax
        mov eax,[ebp+x]
        mov ebx,[ebp+y]
        ;printreg ebx
        ;printreg eax

        movaps xmm1,[eax] ;x[0]
        subps xmm1,[ebx]  ;x[0]-y[0]
        mulps xmm1,xmm1     ;(..)^2
        ;printregps xmm1
        mov edi,[ebp+dddd] ;d
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
        
        mov    edx,[ebp+distance]
        ;printreg edx
        
        movss   [edx],xmm1        ;carico il nuovo valore di distance
        printreg edx
        printregps xmm1
        ;pop     esi
        pop     edi
        pop     ebp
        ret