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

global rowDistance32SdcA

	rowDistance32SdcA:

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
        movups xmm1,[ebx]                   ;c[j*k*sub+i*sub+0]
        subps xmm1,[eax]                    ;c[j*k*sub+i*sub+0] - c[j*k*sub+j_d*sub+0]
        mulps xmm1,xmm1                     ;(..)^2
        ;printregps xmm1
        mov edi,[ebp+subb]                     ;subb
        ;mov ecx,dim
        sub edi,16                          ;subb-16

        mov esi,4                           ;z=4
        xorps xmm2, xmm2
    
    ciclo:
        cmp esi,edi
        jg resto

        movaps xmm0,[ebx+4*esi]             ;c[j*k*sub+i*sub+z]
        ;printregps xmm0
        subps xmm0,[eax+4*esi]              ;c[j*k*sub+i*sub+z] - c[j*k*sub+j_d*sub+z]
        ;printregps xmm0
        mulps xmm0,xmm0                     ;(..)^2
        ;printregps xmm0
        addps xmm1,xmm0                     ;distance+=(..)^2
        ;printregps xmm1
        add esi,4                           ;avanzo di indice

        movaps xmm0,[ebx+4*esi]
        ;printregps xmm0
        subps xmm0,[eax+4*esi]
        ;printregps xmm0
        mulps xmm0,xmm0
        ;printregps xmm0
        addps xmm1,xmm0                     ;distance+=(..)^2
        ;printregps xmm1
        add esi,4

        movaps xmm0,[ebx+4*esi]
        ;printregps xmm0
        subps xmm0,[eax+4*esi]
        ;printregps xmm0
        mulps xmm0,xmm0
        ;printregps xmm0
        addps xmm1,xmm0                     ;distance+=(..)^2
        ;printregps xmm1
        add esi,4

        movaps xmm0,[ebx+4*esi]
        ;printregps xmm0
        subps xmm0,[eax+4*esi]
        ;printregps xmm0
        mulps xmm0,xmm0
        ;printregps xmm0
        addps xmm1,xmm0                     ;distance+=(..)^2
        ;printregps xmm1
        add esi,4

        jmp ciclo
    resto:
        mov edi, [ebp+subb]
        sub edi, 4
    ciclo2:
        cmp esi, edi
        jg resto2
        movaps xmm0,[ebx+4*esi] ;sommo gli ultimi elementi rimanenti
        subps xmm0,[eax+4*esi]
        mulps xmm0,xmm0
        addps xmm2,xmm0
        ;printregps xmm1
        add esi, 4
        jmp ciclo2
    resto2:
        mov edi, [ebp+subb]
    ciclo3:
        cmp esi, edi
        je  fine
        movss xmm0,[ebx+4*esi] ;sommo gli ultimi elementi rimanenti
        subss xmm0,[eax+4*esi]
        mulss xmm0,xmm0
        addss xmm2, xmm0
        ;printregps xmm1
        add esi, 1
        jmp ciclo3
    fine:    
        haddps xmm1,xmm2        ;merge di tutte le somme
        haddps xmm1,xmm1        ;|
        haddps xmm1,xmm1        ;|

        mov eax,[ebp+distance]
        movss  [eax],xmm1

        stop

global rowDistance32SdcU
    rowDistance32SdcU:

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
        ;printregps xmm1
        movups xmm1,[ebx]                   ;c[j*k*sub+i*sub+0]
        ;printregps xmm1
        movups xmm5, [eax]
        subps xmm1,xmm5                    ;c[j*k*sub+i*sub+0] - c[j*k*sub+j_d*sub+0]
        ;printregps xmm1
        mulps xmm1,xmm1                     ;(..)^2
        ;printregps xmm1
        ;mov edi,[ebp+subb]                     ;subb
        ;mov ecx,dim
        sub edi,16                          ;subb-16
        mov esi,4                           ;z=4
        xorps xmm2, xmm2
    
    cicloU:
        cmp esi,edi
        jg restoU

        movups xmm0,[ebx+4*esi]             ;c[j*k*sub+i*sub+z]
        ;printregps xmm0
        movups xmm5, [eax+4*esi]
        subps xmm0,xmm5              ;c[j*k*sub+i*sub+z] - c[j*k*sub+j_d*sub+z]
        ;printregps xmm0
        mulps xmm0,xmm0                     ;(..)^2
        ;printregps xmm0
        addps xmm1,xmm0                     ;distance+=(..)^2
        ;printregps xmm1
        add esi,4                           ;avanzo di indice

        movups xmm0,[ebx+4*esi]
        ;printregps xmm0
        movups xmm5, [eax+4*esi]
        subps xmm0,xmm5
        ;printregps xmm0
        mulps xmm0,xmm0
        ;printregps xmm0
        addps xmm1,xmm0                     ;distance+=(..)^2
        ;printregps xmm1
        add esi,4

        movups xmm0,[ebx+4*esi]
        ;printregps xmm0
        movups xmm5, [eax+4*esi]
        subps xmm0,xmm5
        ;printregps xmm0
        mulps xmm0,xmm0
        ;printregps xmm0
        addps xmm1,xmm0                     ;distance+=(..)^2
        ;printregps xmm1
        add esi,4

        movups xmm0,[ebx+4*esi]
        ;printregps xmm0
        movups xmm5, [eax+4*esi]
        subps xmm0,xmm5
        ;printregps xmm0
        mulps xmm0,xmm0
        ;printregps xmm0
        addps xmm1,xmm0                     ;distance+=(..)^2
        ;printregps xmm1
        add esi,4

        jmp cicloU
    restoU:
        mov edi, [ebp+subb]
        sub edi, 4
    ciclo2U:
        cmp esi, edi
        jg resto2U
        movups xmm0,[ebx+4*esi] ;sommo gli ultimi elementi rimanenti
        movups xmm5, [eax+4*esi]
        subps xmm0,xmm5
        mulps xmm0,xmm0
        addps xmm2,xmm0
        ;printregps xmm1
        add esi, 4
        jmp ciclo2U
    resto2U:
        mov edi, [ebp+subb]
    ciclo3U:
        cmp esi, edi
        je  fineU
        movss xmm0,[ebx+4*esi] ;sommo gli ultimi elementi rimanenti
        subss xmm0,[eax+4*esi]
        mulss xmm0,xmm0
        addss xmm2, xmm0
        ;printregps xmm1
        add esi, 1
        jmp ciclo3U
    fineU:    
        haddps xmm1,xmm2        ;merge di tutte le somme
        haddps xmm1,xmm1        ;|
        haddps xmm1,xmm1        ;|

        mov eax,[ebp+distance]
        movss  [eax],xmm1

        stop