%include "sseutils64.nasm"

section .data

dim equ 4

section .bss

section .text

global dist64A

dist64A:
    start

    vmovaps ymm1,[rdi] ;x[0]
    vsubps ymm1,[rsi]  ;x[0]-y[0]
    vmulps ymm1,ymm1     ;(..)^2
    ;printregyps ymm1
    mov r12,rcx         ;d
    ;mov ecx,dim     ;4
    sub r12,32     ;d-32

    vxorps ymm2, ymm2
    mov r10,8       ;i=4
ciclo:
    cmp r10,r12     ;(j>=d-32)?
    jg resto

    vmovaps ymm0,[rdi+4*r10]             ;uj_x[z]
    vmovaps ymm7,[rdi+4*r10+32]             ;uj_x[z]
    vmovaps ymm8,[rdi+4*r10+64]             ;uj_x[z]
    vmovaps ymm9,[rdi+4*r10+96]             ;uj_x[z]

    vsubps ymm0,[rsi+4*r10]              ;uj_x[z]- c[j*k*subb+i*subb+z
    vsubps ymm7,[rsi+4*r10+32]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    vsubps ymm8,[rsi+4*r10+64]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    vsubps ymm9,[rsi+4*r10+96]              ;uj_x[z]- c[j*k*subb+i*subb+z]
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
    vmovaps ymm0,[rdi+4*r10] ;sommo gli ultimi elementi rimanenti
    vsubps ymm0,[rsi+4*r10]
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
    vmovss  [rdx],xmm1
    stop

global dist64U

dist64U:
    start

    vmovups ymm1,[rdi] ;x[0]
    vsubps ymm1,[rsi]  ;x[0]-y[0]
    vmulps ymm1,ymm1     ;(..)^2
    ;printregyps ymm1
    mov r12,rcx         ;d
    ;mov ecx,dim     ;4
    sub r12,32     ;d-32

    vxorps ymm2, ymm2
    mov r10,8       ;i=4
cicloU:
    cmp r10,r12     ;(j>=d-16)?
    jg restoU

    vmovups ymm0,[rdi+4*r10]             ;x[i]
    vmovups ymm7,[rdi+4*r10+32]           
    vmovups ymm8,[rdi+4*r10+64]           
    vmovups ymm9,[rdi+4*r10+96]           

    vmovups ymm10,[rsi+4*r10]             
    vmovups ymm11,[rsi+4*r10+32]           
    vmovups ymm12,[rsi+4*r10+64]           
    vmovups ymm13,[rsi+4*r10+96]           

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
    vmovups ymm0,[rdi+4*r10] ;sommo gli ultimi elementi rimanenti
    vsubps ymm0,[rsi+4*r10]
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
    vmovss xmm7,[rdi+4*r10] ;sommo gli ultimi elementi rimanenti
    vsubss xmm7,[rsi+4*r10]
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
    vmovss  [rdx],xmm1
    stop

