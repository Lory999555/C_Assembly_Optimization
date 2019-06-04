%include "sseutils.nasm"

section .data
	

    cent	equ		8
    xx      equ     12
    k       equ     16
    ddd     equ     20
    tmp     equ     24  ;risultato di dist
    park    equ     28  ;risultato di centX
    dis     equ     32  ;dis temporanea
    dim1    equ     4

    section .bss
    section .text

    global cent_XA

        cent_XA:
            start
            mov     eax, [ebp+xx]   ;x            
            mov     ebx, [ebp+cent] ;cent
            mov     ecx, [ebp+ddd]  ;d
            mov     edx, [ebp+dis]  ;dis
            mov     edi,    0
            mov     esi, [ebp+k]
            sub     esi, 4
        

            push    eax
            push    ebx
            push    edx
            push    ecx
            push    esi
            push    edi
            call    dist32A
            pop     edi
            pop     esi
            pop     ecx
            pop     edx
            pop     ebx
            pop     eax
            
            xorps   xmm3, xmm3
            movss   xmm3, [edx]
            mov     ecx, [ebp+tmp]      ;tmp
            

        foriA: 
            cmp     edi, esi            ; i> k-4
            jg      fine_A
            
            mov     ebx, [ebp+ddd]
            imul    ebx, edi            ;i*d
            imul    ebx, dim1           ;4*i*d
            add     ebx, [ebp+cent]     ;cent + 4*i*d
            mov     edx, [ebp+ddd]
            push    eax
            push    ebx
            push    ecx
            push    edx
            push    esi
            push    edi
            call    dist32A
            pop     edi
            pop     esi
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax

            ;xorps   xmm4, xmm4      
            movss   xmm4, [ecx]
            comiss  xmm3, xmm4
            jb      avanti1A
            movss   xmm3, xmm4
            mov     ebx, [ebp+park]
            mov     [ebx], edi
        avanti1A:
            add     edi, 1

            mov     ebx, edx
            imul    ebx, edi            ;i*d
            imul    ebx, dim1           ;4*i*d
            add     ebx, [ebp+cent]     ;cent + 4*i*d
            push    eax
            push    ebx
            push    ecx
            push    edx
            push    esi
            push    edi
            call    dist32A
            pop     edi
            pop     esi
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax

            ;xorps   xmm4, xmm4          ;;;;;;;;;;;;;;;;;;;;;;;;
            movss   xmm4, [ecx]
            comiss  xmm3, xmm4
            jb      avanti2A
            movss   xmm3, xmm4
            mov     ebx, [ebp+park]
            mov     [ebx], edi
        avanti2A:
            add     edi, 1

            mov     ebx, edx
            imul    ebx, edi            ;i*d
            imul    ebx, dim1           ;4*i*d
            add     ebx, [ebp+cent]     ;cent + 4*i*d
            push    eax
            push    ebx
            push    ecx
            push    edx
            push    esi
            push    edi
            call    dist32A
            pop     edi
            pop     esi
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax

            ;xorps   xmm4, xmm4          ;;;;;;;;;;;;;;;;;;;;;;;;
            movss   xmm4, [ecx]
            comiss  xmm3, xmm4
            jb      avanti3A
            movss   xmm3, xmm4
            mov     ebx, [ebp+park]
            mov     [ebx], edi

        avanti3A:
            add     edi, 1

            mov     ebx, edx
            imul    ebx, edi            ;i*d
            imul    ebx, dim1           ;4*i*d
            add     ebx, [ebp+cent]     ;cent + 4*i*d
            push    eax
            push    ebx
            push    ecx
            push    edx
            push    esi
            push    edi
            call    dist32A
            pop     edi
            pop     esi
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax

            ;xorps   xmm4, xmm4          ;;;;;;;;;;;;;;;;;;;;;;;;
            movss   xmm4, [ecx]
            comiss  xmm3, xmm4
            jb      avanti4A
            movss   xmm3, xmm4
            mov     ebx, [ebp+park]
            mov     [ebx], edi
        avanti4A:
            add     edi, 1         
            jmp     foriA
        fine_A:
            stop


    global cent_XU

        cent_XU:
            start
            mov     eax, [ebp+xx]   ;x            
            mov     ebx, [ebp+cent] ;cent
            mov     ecx, [ebp+ddd]  ;d
            mov     edx, [ebp+dis]  ;dis
            mov     edi,    0
            mov     esi, [ebp+k]
            sub     esi, 4
        

            push    eax
            push    ebx
            push    edx
            push    ecx
            push    esi
            push    edi
            call    dist32U
            pop     edi
            pop     esi
            pop     ecx
            pop     edx
            pop     ebx
            pop     eax
            
            xorps   xmm3, xmm3
            movss   xmm3, [edx]
            mov     ecx, [ebp+tmp]      ;tmp
            

        foriU: 

            cmp     edi, esi            ; i> k-4
            jg      fine_U
            
            mov     ebx, [ebp+ddd]
            imul    ebx, edi            ;i*d
            imul    ebx, dim1           ;4*i*d
            add     ebx, [ebp+cent]     ;cent + 4*i*d
            mov     edx, [ebp+ddd]
            push    eax
            push    ebx
            push    ecx
            push    edx
            push    esi
            push    edi
            call    dist32U
            pop     edi
            pop     esi
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax

            ;xorps   xmm4, xmm4      
            movss   xmm4, [ecx]
            comiss  xmm3, xmm4
            jb      avanti1U
            movss   xmm3, xmm4
            mov     ebx, [ebp+park]
            mov     [ebx], edi
        avanti1U:
            add     edi, 1

            mov     ebx, edx
            imul    ebx, edi            ;i*d
            imul    ebx, dim1           ;4*i*d
            add     ebx, [ebp+cent]     ;cent + 4*i*d
            push    eax
            push    ebx
            push    ecx
            push    edx
            push    esi
            push    edi
            call    dist32U
            pop     edi
            pop     esi
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax

            ;xorps   xmm4, xmm4          ;;;;;;;;;;;;;;;;;;;;;;;;
            movss   xmm4, [ecx]
            comiss  xmm3, xmm4
            jb      avanti2U
            movss   xmm3, xmm4
            mov     ebx, [ebp+park]
            mov     [ebx], edi
        avanti2U:
            add     edi, 1

            mov     ebx, edx
            imul    ebx, edi            ;i*d
            imul    ebx, dim1           ;4*i*d
            add     ebx, [ebp+cent]     ;cent + 4*i*d
            push    eax
            push    ebx
            push    ecx
            push    edx
            push    esi
            push    edi
            call    dist32U
            pop     edi
            pop     esi
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax

            ;xorps   xmm4, xmm4          ;;;;;;;;;;;;;;;;;;;;;;;;
            movss   xmm4, [ecx]
            comiss  xmm3, xmm4
            jb      avanti3U
            movss   xmm3, xmm4
            mov     ebx, [ebp+park]
            mov     [ebx], edi

        avanti3U:
            add     edi, 1

            mov     ebx, edx
            imul    ebx, edi            ;i*d
            imul    ebx, dim1           ;4*i*d
            add     ebx, [ebp+cent]     ;cent + 4*i*d
            push    eax
            push    ebx
            push    ecx
            push    edx
            push    esi
            push    edi
            call    dist32U
            pop     edi
            pop     esi
            pop     edx
            pop     ecx
            pop     ebx
            pop     eax

            ;xorps   xmm4, xmm4          ;;;;;;;;;;;;;;;;;;;;;;;;
            movss   xmm4, [ecx]
            comiss  xmm3, xmm4
            jb      avanti4U
            movss   xmm3, xmm4
            mov     ebx, [ebp+park]
            mov     [ebx], edi
        avanti4U:
            add     edi, 1         
            jmp     foriU
        fine_U:
            stop




section .data
    x		    equ		28
    y           equ     24
    distance	equ		20
    dddd        equ		16
    dim         equ     4

    section .bss
    section .text

    dist32U:
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
        mov edi,[ebp+dddd] ;d
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
        mov edi, [ebp+dddd]
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
        mov edi, [ebp+dddd]
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







dist32A:
        push ebp
        mov ebp, esp

        mov eax,[ebp+x]
        mov ebx,[ebp+y]

        movaps xmm1,[eax] ;x[0]
        subps xmm1,[ebx]  ;x[0]-y[0]
        mulps xmm1,xmm1     ;(..)^2
        ;printregps xmm1
        mov edi,[ebp+dddd] ;d
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
        mov edi, [ebp+dddd]
        sub edi, 4
    ciclo2:
        cmp esi, edi
        jg  resto2
        movaps xmm0,[eax+4*esi] ;sommo gli ultimi elementi rimanenti
        subps xmm0,[ebx+4*esi]
        mulps xmm0,xmm0
        addps xmm2, xmm0
        add esi, 4
        jmp ciclo2

    resto2:
        mov edi, [ebp+dddd]
    ciclo3:
        cmp esi, edi
        je  fineA
        movss xmm0,[eax+4*esi] ;sommo gli ultimi elementi rimanenti
        subss xmm0,[ebx+4*esi]
        mulss xmm0,xmm0
        addss xmm2, xmm0
        add esi, 1
        jmp ciclo3
    fineA:
        haddps xmm1,xmm2        ;merge di tutte le somme
        haddps xmm1,xmm1        ;|
        haddps xmm1,xmm1        ;|
        ;printregps xmm1        

        mov eax,[ebp+distance]
        movss [eax],xmm1        ;carico il nuovo valore di distance
        pop     ebp
        ret