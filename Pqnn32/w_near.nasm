%include "sseutils.nasm"

section .data

    x1          equ     8
    cent        equ     12
    tmp     	equ 	16
    d1		    equ     20
    res_w       equ     24
    res_d       equ     28
    w           equ     32
    n           equ     36
    max         equ     40
    dim1        equ     4

    extern max_heap
    section .bss   
    section .text
 
    global w_near
        w_near:
            start
            ;printreg eax
            mov eax, [ebp+x1]
            mov ecx, [ebp+tmp]
            mov edx, [ebp+d1]
            mov esi, [ebp+w]
            mov edi, 0
            xorps xmm5, xmm5

        fori:
            cmp     edi, esi        ; i> w
            jge      avanti2
            ;printreg edi
            mov     ebx, [ebp+d1]
            imul    ebx, edi            ;i*d
            imul    ebx, dim1           ;4*i*d
            add     ebx, [ebp+cent]     ;cent + 4*i*d
            push    eax
            push    ebx
            push    ecx
            push    edx
            push    edi
            push    esi
            call    dist32
            pop     esi
            pop     edi
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax

            mov     ebx, edi
            imul    ebx, dim1
            add     ebx, [ebp+res_w]
            mov     [ebx], edi          ;reult_w[i]=i
            ;printreg ebx
            ;printreg edi
            ;printreg [ebx]
            mov     ebx, edi
            imul    ebx, dim1
            add     ebx, [ebp+res_d]
            movss   xmm4, [ecx]
            movss   [ebx], xmm4         ;reult_d[i]=tmp
            
            comiss  xmm4, xmm5
            jb      avanti1
            movss   xmm5, xmm4
            ;printregps xmm5

        avanti1:
            add     edi, 1
            jmp     fori
        
        avanti2:
            mov     edi,[ebp+w]
            mov     esi, [ebp+n]
        forj:
            cmp     edi, esi        ; i>n
            jge      forj
            mov     ebx, [ebp+d1]
            imul    ebx, edi            ;i*d
            imul    ebx, dim1           ;4*i*d
            add     ebx, [ebp+cent]     ;cent + 4*i*d
            push    eax
            push    ebx
            push    ecx
            push    edx
            push    edi
            push    esi
            call    dist32
            pop     esi
            pop     edi
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax
            
            movss   xmm4, [ecx]
            printregps xmm4
            printregps xmm5
            comiss  xmm5, xmm4
            jb      avanti3  
            mov     eax, [ebp+max]   
            push    eax
            push    1 
            mov     ebx, [ebp+w]
            push    ebx
            extractps ebx, xmm5, 0
            push    ebx
            mov     ebx, [ecx]
            push    ebx
            push    edi
            mov     ebx, [ebp+res_d]
            push    ebx
            mov     ebx, [ebp+res_w]
            push    ebx
            call    max_heap
            pop     ebx
            pop     ebx
            pop     edi
            pop     ecx
            pop     ebx
            pop     ebx
            pop     ebx
            pop     eax
            movss   xmm5,[eax]
            printreg eax
            printregps xmm5
            
            ;printreg [eax]
            
        avanti3:
            mov     eax, [ebp+x1]
            mov     ecx, [ebp+tmp]
            add     edi, 1
            jmp     forj

        
        fine:
            stop

                









section .data
    x		        equ	28
    y                   equ     24
    distance	        equ	20
    d2                  equ	16
    dim                 equ     4

    section .bss
    section .text

    dist32:
        ;printreg [ebp+y]
        ;printreg [ebp+distance]
        ;printreg [ebp+dddd]


        push ebp
        mov ebp, esp
        mov eax,[ebp+x]
        mov ebx,[ebp+y]
        xorps xmm1, xmm1
        movaps xmm1,[eax] ;x[0]
        movaps xmm7, [ebx]
        subps xmm1,xmm7  ;x[0]-y[0]
        mulps xmm1,xmm1     ;(..)^2
        ;printregps xmm1
        mov edi,[ebp+d2] ;d
        ;mov ecx,dim     ;4
        sub edi,16     ;d-16
        xorps xmm2, xmm2
        mov esi,4       ;i=4
    cicloU:
        cmp esi,edi     ;(j>=d-16)?
        jg restoU
        movaps xmm0,[eax+4*esi] ;x[i]
        movaps xmm7, [ebx+4*esi]
        subps xmm0,xmm7  ;x[i]-y[i]
        mulps xmm0,xmm0         ;(..)^2
        ;printregps xmm0edxedx
        addps xmm1,xmm0         ;distance+=(..)^2
        ;printregps xmm1
        add esi,4               ;avanzo di indice

        movaps xmm0,[eax+4*esi] ;x[i]
        movaps xmm7, [ebx+4*esi]
        subps xmm0,xmm7  ;x[i]-y[i]
        mulps xmm0,xmm0         ;(..)^2
        ;printregps xmm0
        addps xmm1,xmm0         ;distance+=(..)^2
        ;printregps xmm1
        add esi,4               ;avanzo di indice

        movaps xmm0,[eax+4*esi] ;x[i]
        movaps xmm7, [ebx+4*esi]
        subps xmm0,xmm7  ;x[i]-y[i]
        mulps xmm0,xmm0         ;(..)^2
        ;printregps xmm0
        addps xmm1,xmm0         ;distance+=(..)^2
        ;printregps xmm1
        add esi,4               ;avanzo di indice

        movaps xmm0,[eax+4*esi] ;x[i]
        movaps xmm7, [ebx+4*esi]
        subps xmm0,xmm7  ;x[i]-y[i]
        mulps xmm0,xmm0         ;(..)^2
        ;printregps xmm0
        addps xmm1,xmm0         ;distance+=(..)^2
        ;printregps xmm1
        add esi,4               ;avanzo di indice

        jmp cicloU
    restoU:
        mov edi, [ebp+d2]
        sub edi, 4
    ciclo2U:
        cmp esi, edi
        jg  resto2U
        movaps xmm0,[eax+4*esi] ;sommo gli ultimi elementi rimanenti
        movaps xmm7, [ebx+4*esi]
        subps xmm0,xmm7
        mulps xmm0,xmm0
        addps xmm2, xmm0
        add esi, 4
        jmp ciclo2U

    resto2U:
        mov edi, [ebp+d2]
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

        mov eax,[ebp+distance]
        movss [eax],xmm1        ;carico il nuovo valore di distance        pop     edi
        pop     ebp
        ret