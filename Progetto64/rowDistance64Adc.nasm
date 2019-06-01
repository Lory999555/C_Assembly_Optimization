%include "sseutils64.nasm"

section .data

subb equ 16

dim equ 4


section .bss

section .text

global rowDistance64AdcA

global rowDistance64AdcU

rowDistance64AdcA:

    push		rbp				; salva il Base Pointer
    mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
    pushaq						; salva i registri generali

    mov r11,[rbp+subb]
    ;vprintreg r11

    imul r8,dim                         ;prendo j
    imul r8,r11                        ;j*sub
    imul r8,rcx                        ;j*k*subb
    ;printreg edx

                            
    imul rcx,dim                        ;prendo i
    imul rcx,r11                        ;i*subb
    add rcx,r8                          ;j*k*subb+i*subb
    add rcx,rdx                         ;c[j*k*subb+i*subb]

    ;vxorps xmm5,xmm5                     ;distance=0 

    vmovaps ymm1,[rsi]                   ;uj_x[0]
    vsubps ymm1,[rcx]                    ;uj_x[0]-c[0*k*subb+i*subb]
    vmulps ymm1,ymm1                     ;(..)^2
    printregyps ymm1
    ;mov r11,[rbp+subb]                     ;subb
    sub r11,32                          ;subb-16

    mov r12,8                            ;z=8
    vxorps ymm2,ymm2
cicloAAA:
    cmp r12,r11
    jge fine

    vmovaps ymm0,[rsi+4*r12]             ;uj_x[z]
    printregyps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice

    vmovaps ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice

    vmovaps ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice

    vmovaps ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice

    jmp cicloAAA
resto:
    mov r11,[rbp+subb]                 ;subb
    sub r11,4                          ;subb-4
ciclo2:
    cmp r12, r11
    jg resto2
    vmovaps ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm2,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice
    jmp ciclo2
resto2:
    mov r11, [rbp+subb]
ciclo3:
    cmp r12, r11
    je  fine
    movss xmm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    subss xmm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    mulss xmm0,xmm0                     ;(..)^2
    ;printregps ymm0
    addss xmm2,xmm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,1                           ;avanzo di indice
    jmp ciclo3
fine:
    printregyps ymm1
    printregyps ymm2
    vhaddps ymm1,ymm2        ;merge di tutte le somme
    vhaddps ymm1,ymm1        ;|
    vhaddps ymm1,ymm1        ;|

    movss  [rdx],xmm1

    stop

rowDistance64AdcU:

    push		rbp				; salva il Base Pointer
    mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
    pushaq						; salva i registri generali

    mov r11,[rbp+subb]


    imul r8,dim                         ;prendo j
    imul r8,r11                        ;j*sub
    imul r8,rcx                        ;j*k*subb
    ;printreg edx

                            
    imul rcx,dim                        ;prendo i
    imul rcx,r11                        ;i*subb
    add rcx,r8                         ;j*k*subb+i*subb
    add rcx,rdx                         ;c[j*k*subb+i*subb]

    ;vxorps xmm5,xmm5                     ;distance=0 

    vmovups ymm1,[rsi]                   ;uj_x[0]
    vsubps ymm1,[rcx]                    ;uj_x[0]-c[0*k*subb+i*subb]
    vmulps ymm1,ymm1                     ;(..)^2
    ;printregps ymm1
    ;mov r11,[rbp+subb]                     ;subb
    sub r11,16                          ;subb-16

    mov r12,8                            ;z=8
    vxorps ymm2,ymm2
cicloU:
    cmp r12,r11
    jge fineU

    vmovups ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice

    vmovups ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice

    vmovups ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice

    vmovups ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice

    jmp cicloU
restoU:
    mov r11,[rbp+subb]                 ;subb
    sub r11,8                          ;subb-8
ciclo2U:
    cmp r12, r11
    jg resto2U
    vmovups ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm2,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice
    jmp ciclo2U
resto2U:
    mov r11, [rbp+subb]
ciclo3U:
    cmp r12, r11
    je  fineU
    movss xmm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    subss xmm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    mulss xmm0,xmm0                     ;(..)^2
    ;printregps ymm0
    addss xmm2,xmm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,1                           ;avanzo di indice
    jmp ciclo3U
fineU:
    
    vhaddps ymm1,ymm2        ;merge di tutte le somme
    vhaddps ymm1,ymm1        ;|
    vhaddps ymm1,ymm1        ;|

    movss  [rdx],xmm1

    stop