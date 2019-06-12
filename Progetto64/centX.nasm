%include "sseutils64.nasm"

section .data

dis     equ     16  ;dis temporanea
dim1    equ     4

section .bss
section .text

global cent_XA

cent_XA:
    start
    mov     r14, [rbp+dis]  ;dis

    mov     r11, 0
    sub     rdx, 4      ;k-4
    mov     rax,rdi

    ;pushaq
    push    rax
    push    r13
    push    r12
    push    r10
    call    dist64A
    pop     r10
    pop     r12
    pop     r13
    pop     rax
    ;popaq

    
    vxorps   xmm3, xmm3
    vmovss   xmm3, [r14]         ;dis
foriA: 
    cmp     r11, rdx            ; i> k-4
    jg      fine_A
    ;vprintreg r11

    mov     r13, rcx            ;d
    imul    r13,r11             ;i*d
    imul    r13,dim1            ;4*i*d
    ;vprintreg rax
    ;vprintreg r13
    mov     rax,rdi
    add     rax,r13            ;cent[4*i*d]
    ;vprintreg rax
    ;pushaq
    push    rax
    push    r13
    push    r12
    push    r10
    call    dist64A
    pop     r10
    pop     r12
    pop     r13
    pop     rax
    ;popaq

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti1A
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i

avanti1A:
    add     r11, 1

    ;vprintreg r13
    ;vprintreg rax
    mov     r13, rcx            ;d
    imul    r13,r11             ;4*i*d
    imul    r13,dim1
    mov     rax,rdi
    add     rax,r13            ;cent[4*i*d]
    ;vprintreg rax

    ;pushaq
    push    rax
    push    r13
    push    r12
    push    r10
    call    dist64A
    pop     r10
    pop     r12
    pop     r13
    pop     rax
    ;popaq

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti2A
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i
avanti2A:
    add     r11, 1

    mov     r13,rcx            ;d
    imul    r13,r11             ;4*i*d
    imul    r13,dim1
    mov     rax,rdi
    add     rax,r13            ;cent[4*i*d]
    ;vprintreg rax

    ;pushaq
    push    rax
    push    r13
    push    r12
    push    r10
    call    dist64A
    pop     r10
    pop     r12
    pop     r13
    pop     rax
    ;popaq

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti3A
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i

avanti3A:
    add     r11, 1

    mov     r13, rcx            ;d
    imul    r13,r11             ;4*i*d
    imul    r13,dim1
    mov     rax,rdi
    add     rax,r13            ;cent[4*i*d]
   ;vprintreg rax

    ;pushaq
    push    rax
    push    r13
    push    r12
    push    r10
    call    dist64A
    pop     r10
    pop     r12
    pop     r13
    pop     rax
    ;popaq

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti4A
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i
avanti4A:
    add     r11, 1         
    jmp     foriA
fine_A:
    stop


global cent_XU

cent_XU:
    start
    mov     r14, [rbp+dis]  ;dis

    mov     r11, 0
    sub     rdx, 4      ;k-4
    mov     rax,rdi

    ;pushaq
    push    rax
    push    r13
    push    r12
    push    r10
    call    dist64U
    pop     r10
    pop     r12
    pop     r13
    pop     rax
    ;popaq

    
    vxorps   xmm3, xmm3
    vmovss   xmm3, [r14]         ;dis
foriU: 
    cmp     r11, rdx            ; i> k-4
    jg      fine_U
    ;vprintreg r11

    mov     r13, rcx            ;d
    imul    r13,r11             ;i*d
    imul    r13,dim1            ;4*i*d
    ;vprintreg rax
    ;vprintreg r13
    mov     rax,rdi
    add     rax,r13            ;cent[4*i*d]
    ;vprintreg rax
    ;pushaq
    push    rax
    push    r13
    push    r12
    push    r10
    call    dist64U
    pop     r10
    pop     r12
    pop     r13
    pop     rax
    ;popaq

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti1U
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i

avanti1U:
    add     r11, 1

    ;vprintreg r13
    ;vprintreg rax
    mov     r13, rcx            ;d
    imul    r13,r11             ;4*i*d
    imul    r13,dim1
    mov     rax,rdi
    add     rax,r13            ;cent[4*i*d]
    ;vprintreg rax

    ;pushaq
    push    rax
    push    r13
    push    r12
    push    r10
    call    dist64U
    pop     r10
    pop     r12
    pop     r13
    pop     rax
    ;popaq

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti2U
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i
avanti2U:
    add     r11, 1

    mov     r13,rcx            ;d
    imul    r13,r11             ;4*i*d
    imul    r13,dim1
    mov     rax,rdi
    add     rax,r13            ;cent[4*i*d]
    ;vprintreg rax

    ;pushaq
    push    rax
    push    r13
    push    r12
    push    r10
    call    dist64U
    pop     r10
    pop     r12
    pop     r13
    pop     rax
    ;popaq

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti3U
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i

avanti3U:
    add     r11, 1

    mov     r13, rcx            ;d
    imul    r13,r11             ;4*i*d
    imul    r13,dim1
    mov     rax,rdi
    add     rax,r13            ;cent[4*i*d]
   ;vprintreg rax

    ;pushaq
    push    rax
    push    r13
    push    r12
    push    r10
    call    dist64U
    pop     r10
    pop     r12
    pop     r13
    pop     rax
    ;popaq

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti4U
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i
avanti4U:
    add     r11, 1         
    jmp     foriU
fine_U:
    stop

section .data

dim equ 4

section .bss

section .text

dist64A:
    start

    vmovaps ymm1,[rsi] ;x[0]
    vsubps ymm1,[rax]  ;x[0]-y[0]
    vmulps ymm1,ymm1     ;(..)^2
    ;printregyps ymm1
    mov r12,rcx         ;d
    ;mov rcx,dim     ;4
    sub r12,32     ;d-16

    vxorps ymm2, ymm2
    mov r10,8       ;j=8
ciclo:
    cmp r10,r12     ;(j>=d-16)?
    jg resto

    vmovaps ymm0,[rsi+4*r10]             ;uj_x[z]
    vmovaps ymm7,[rsi+4*r10+32]             ;uj_x[z]
    vmovaps ymm8,[rsi+4*r10+64]             ;uj_x[z]
    vmovaps ymm9,[rsi+4*r10+96]             ;uj_x[z]

    vsubps ymm0,[rax+4*r10]              ;uj_x[z]- c[j*k*subb+i*subb+z
    vsubps ymm7,[rax+4*r10+32]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    vsubps ymm8,[rax+4*r10+64]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    vsubps ymm9,[rax+4*r10+96]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    vmulps ymm7,ymm7                     ;(..)^2
    vmulps ymm8,ymm8                     ;(..)^2
    vmulps ymm9,ymm9                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    vaddps ymm1,ymm7                     ;distance+=(..)^2
    vaddps ymm1,ymm8                     ;distance+=(..)^2
    vaddps ymm1,ymm9                     ;distance+=(..)^2
    ;printregps ymm1
    add r10,32                           ;avanzo di indice

    jmp ciclo
resto:
    mov r12, rcx
    sub r12, 8
ciclo2:
    cmp r10, r12
    jg  fine
    vmovaps ymm0, [rsi+4*r10] ;sommo gli ultimi elementi rimanenti
    vsubps ymm0,[rax+4*r10]
    vmulps ymm0,ymm0
    vaddps ymm2, ymm0
    add r10, 8
    jmp ciclo2
fine:
    vaddps ymm1,ymm2        ;merge di tutte le somme
    ;printregyps ymm1
    vhaddps ymm1,ymm1        ;|
    ;printregyps ymm1
    vhaddps ymm1,ymm1        ;|
    ;printregyps ymm1
    vperm2f128 ymm5,ymm1,ymm1,1
    vaddss xmm1,xmm5
    ;printregyps ymm1
    vmovss  [r14],xmm1
    stop

dist64U:
    start

    vmovups ymm1,[rsi] ;x[0]
    vsubps ymm1,[rax]  ;x[0]-y[0]
    vmulps ymm1,ymm1     ;(..)^2
    ;printregyps ymm1
    mov r12,rcx         ;d
    ;mov rcx,dim     ;4
    sub r12,32     ;d-16

    vxorps ymm2, ymm2
    mov r10,8       ;i=4
cicloU:
    cmp r10,r12     ;(j>=d-16)?
    jg restoU

    vmovups ymm0,[rsi+4*r10]             ;x[i]
    vmovups ymm7,[rsi+4*r10+32]           
    vmovups ymm8,[rsi+4*r10+64]           
    vmovups ymm9,[rsi+4*r10+96]           

    vmovups ymm10,[rax+4*r10]             
    vmovups ymm11,[rax+4*r10+32]           
    vmovups ymm12,[rax+4*r10+64]           
    vmovups ymm13,[rax+4*r10+96]           

    vsubps ymm0,ymm10                   ;x[i]-y[i]
    vsubps ymm7,ymm11              
    vsubps ymm8,ymm12              
    vsubps ymm9,ymm13              
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    vmulps ymm7,ymm7                     ;(..)^2
    vmulps ymm8,ymm8                     ;(..)^2
    vmulps ymm9,ymm9                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    vaddps ymm1,ymm7                     ;distance+=(..)^2
    vaddps ymm1,ymm8                     ;distance+=(..)^2
    vaddps ymm1,ymm9                     ;distance+=(..)^2
    ;printregps ymm1
    add r10,32                           ;avanzo di indice

    jmp cicloU
restoU:
    mov r12, rcx
    sub r12, 8
ciclo2U:
    cmp r10, r12
    jg  resto2U
    vmovups ymm0, [rsi+4*r10] ;sommo gli ultimi elementi rimanenti
    vsubps ymm0,[rax+4*r10]
    vmulps ymm0,ymm0
    vaddps ymm2, ymm0
    ;printregyps ymm2
    add r10, 8
    jmp ciclo2U

resto2U:
    mov r12, rcx
ciclo3U:
    cmp r10, r12
    je  fineU
    vmovss xmm7, [rsi+4*r10] ;sommo gli ultimi elementi rimanenti
    vsubss xmm7,[rax+4*r10]
    vmulss xmm7, xmm7
    ;printregyps ymm7
    vaddps ymm2, ymm7
    ;printregyps ymm2
    add r10, 1
    jmp ciclo3U
fineU:
    vaddps ymm1,ymm2        ;merge di tutte le somme
    ;printregyps ymm1
    vhaddps ymm1,ymm1        ;|
    ;printregyps ymm1
    vhaddps ymm1,ymm1        ;|
    ;printregyps ymm1
    vperm2f128 ymm5,ymm1,ymm1,1
    vaddss xmm1,xmm5
    ;printregyps ymm1
    vmovss  [r14],xmm1
    stop